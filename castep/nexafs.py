#!/usr/bin/env python3

import os
from functools import lru_cache, partial
from multiprocessing import Pool
from typing import List, Tuple

import click
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


class Nexafs:
    def __init__(
        self,
        n_procs,
        begin,
        end,
        initial_peak,
        omega_1,
        omega_2,
        mix_1,
        mix_2,
        molecule,
        surface,
        excited_atom,
        index_start,
        index_end,
        spectrum_type,
        phi,
        theta,
        excited_nth_element,
        plot_i_atoms,
    ):
        """
        Parameters
        ----------
            n_procs : int
                number of processors to use
            begin : float
                start of energy range to plot
            end : float
                end of energy range to plot
            initial_peak : float
                energy value of the initial peak
            omega_1 : float
                full width at half maximum for region 1
            omega_2 : float
                full width at half maximum for region 2
            mix_1 : float
                Gaussian-Lorentzian mixing parameter for region 1
            mix_2 : float
                Gaussian-Lorentzian mixing parameter for region 2
            molecule : str
                adsorbate name
            surface : str
                surface element
            excited_atom : str
                element symbol of the excited atom
            index_start : int
                index of the first excited atom
            index_end : int
                index of the last excited atom
            spectrum_type : str
                type of spectrum to plot
            phi : list
                phi angles
            theta : list
                theta angles
            excited_nth_element : int
                the element index of the excited atom of all the elements in the system
                - this will always be the last element in the system, so if system
                contains C, Cu, C:exc it will be 3
            plot_i_atoms : bool
                plot individual atom spectra
        """
        # Number of processors to use
        self.n_procs = n_procs

        # Broadening parameters
        # Start and end values of spectra
        self.start = begin
        self.stop = end

        # Energy of the lead peak
        self.first_peak = initial_peak
        self.omega_1 = omega_1
        self.omega_2 = omega_2
        self.mix_1 = mix_1
        self.mix_2 = mix_2

        # Details about the chemical system
        self.molecule = molecule
        self.metal = surface
        self.element = excited_atom

        # Index range of the atom directories created by autoscript.py
        # self.start_atom = index_start
        # self.end_atom = index_end
        start_atom = index_start
        end_atom = index_end

        # Type of NEXAFS spectrum to output
        # 1 for total summed NEXAFS, 2 for angular, 3 for polarised, and 4 for
        # average polarised
        spectrum_type_opts = np.array(
            ["tot_summed", "angular", "polarised", "avg_polarised"]
        )
        self.n_type = np.where(spectrum_type_opts == spectrum_type)[0][0] + 1

        # Convert phi and theta to numpy arrays
        self.phi = np.array(phi)
        self.theta = np.array(theta)

        self.atom = excited_nth_element
        self.plot_i_atoms = plot_i_atoms

        # Set the start and end point of the linearly increasing broadening
        # change from the first and last ranges with respect to the leading
        # peak
        self.ewid_1 = initial_peak + 5.0
        self.ewid_2 = initial_peak + 15.0

        # Create a list of all the atoms
        atom_ids = list(range(start_atom, end_atom + 1))

        # Set up a list of all the directories all the data is in C48/, C49/... C57/
        self.dirs = np.array([])
        for n in atom_ids:
            direc = f"{self.element}{str(n)}/"
            if os.path.isdir(direc):
                self.dirs = np.append(self.dirs, direc)

        # Deltas file for each molpdos run
        self.fname = f"/{self.molecule}_{self.metal}_{self.atom}_1_1_1_deltas.dat"

        # Get the length of the deltas file
        tmp_bands = np.loadtxt(
            f"{self.element}{str(atom_ids[0])}/t{self.theta[0]}_p{self.phi[0]}{self.fname}"
        )
        tmp_bands_l = len(tmp_bands)

        # Create arrays with sizes of the system to use
        self.peaks = np.zeros([len(atom_ids), tmp_bands_l])
        self.bands = np.zeros([len(atom_ids), tmp_bands_l])

    def get_nexafs_data(
        self,
        get_i_atom=False,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Parse the NEXAFS data for all given theta and phi angles.

        Parameters
        ----------
            peaks : np.ndarray
                array to store the Dirac peak positions
            bands : np.ndarray
                array to store the Dirac peak intensities
            get_i_atom : bool
                whether to parse and save the individual atom spectra to a file

        Returns
        -------
            flat_peaks : np.ndarray
                flattened array of Dirac peak positions
            flat_bands : np.ndarray
                flattened array of Dirac peak intensities
            peaks : np.ndarray
                array of Dirac peak positions
            bands : np.ndarray
                array of Dirac peak intensities
        """

        print(f"Parsing NEXAFS data for self.theta={self.theta} self.phi={self.phi}...")
        for i in range(len(self.dirs)):
            nexafs = np.loadtxt(f"{self.dirs[i]}t{self.theta}_p{self.phi}{self.fname}")
            self.peaks[i, :] = nexafs[:, 0]
            self.bands[i, :] = nexafs[:, self.n_type]

            if get_i_atom:  # Save each atom spectrum to a file
                i_atom_peaks = nexafs[:, 0]
                i_atom_bands = nexafs[:, self.n_type]
                np.savetxt(
                    f"{self.dirs[i]}t{self.theta}_p{self.phi}/{self.molecule}_"
                    f"{self.metal}_deltas_t{self.theta}_p{self.phi}.txt",
                    np.vstack((i_atom_peaks, i_atom_bands)).T,
                    header="# <x in eV> Intensity",
                )

        flat_peaks = self.peaks.flatten()
        flat_bands = self.bands.flatten()

        # Write out all of the data into a delta peaks file
        np.savetxt(
            f"{self.molecule}_{self.metal}_deltas_t{self.theta}_p{self.phi}.txt",
            np.vstack((flat_peaks, flat_bands)).T,
            header="# <x in eV> Intensity",
        )

        print(
            f"Finished writing out delta peaks file for theta={self.theta} phi={self.phi}"
        )

        return flat_peaks, flat_bands, self.peaks, self.bands

    @staticmethod
    @lru_cache(maxsize=128)
    def _schmid_pseudo_voigt(domain, m, E, omega) -> np.ndarray:
        """
        Apply broadening scheme for XPS spectra. See the following reference for more details:
        https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/sia.5521

        Parameters
        ----------
            domain : np.ndarray
                range of x-values to include in the spectrum
            m : float
                Gaussian-Lorentzian mixing parameter
            E : float
                line centre/dirac peak
            omega : float
                full width at half maximum (omega > 0)

        Returns
        -------
            data : np.ndarray
                broadened spectrum
        """

        return (1 - m) * np.sqrt((4 * np.log(2)) / (np.pi * omega**2)) * np.exp(
            -(4 * np.log(2) / omega**2) * (domain - E) ** 2
        ) + m * (1 / (2 * np.pi)) * (omega / ((omega / 2) ** 2 + (domain - E) ** 2))

    @staticmethod
    def _normalise(data, domain, ewid_1, norm_val=None) -> tuple[np.ndarray, float]:
        """
        Normalise spectrum such that first peak has intensity of 1.

        Parameters
        ----------
            data : np.ndarray
                broadened spectrum to normalise
            domain : np.ndarray
                range of x-values to include in the spectrum
            ewid_1 : float
                x-value at which the first peak region ends
            norm_val : float
                value to normalise the spectrum by

        Returns
        -------
            data : np.ndarray
                normalised spectrum
            k_edge_max : float
                maximum intensity of the first peak
        """
        # TODO: Change norm_val/k_edge_max to use getter/setter methods

        # Find height of first peak
        if norm_val is None:
            k_edge_max = np.max(data[np.where(domain <= ewid_1)])
        else:
            k_edge_max = norm_val

        # Scale the rest of the data proportionally
        data /= k_edge_max

        return data, k_edge_max

    @staticmethod
    def _mp_broaden(domain, mixing, dirac_peaks, sigma, coeffs) -> np.ndarray:
        """
        Perform the broadening in parallel.

        Parameters
        ----------
            domain : np.ndarray
                range of x-values to include in the spectrum
            mixing : np.ndarray
                Gaussian-Lorentzian mixing parameter
            dirac_peaks : np.ndarray
                line centre/dirac peak
            sigma : np.ndarray
                linearly interpolated full width at half maximum
            coeffs : np.ndarray
                Dirac peak intensities

        Returns
        -------
            data : np.ndarray
                broadened spectrum
        """
        data = np.sum(
            Nexafs._schmid_pseudo_voigt(domain, mixing, dirac_peaks, sigma) * coeffs
        )

        return data

    @staticmethod
    def _sp_broaden(data, domain, mixing, dirac_peaks, sigma, coeffs) -> np.ndarray:
        """
        Perform the broadening in serial.

        Parameters
        ----------
            data : np.ndarray
                non-broadened spectrum
            domain : np.ndarray
                range of x-values to include in the spectrum
            mixing : np.ndarray
                Gaussian-Lorentzian mixing parameter
            dirac_peaks : np.ndarray
                line centre/dirac peak
            sigma : np.ndarray
                linearly interpolated full width at half maximum
            coeffs : np.ndarray
                Dirac peak intensities

        Returns
        -------
            data : np.ndarray
                broadened spectrum
        """
        for i in tqdm(range(len(domain))):
            data[i] = np.sum(
                Nexafs._schmid_pseudo_voigt(domain[i], mixing, dirac_peaks, sigma)
                * coeffs
            )

        return data

    def broaden(
        self,
        dirac_peaks,
        coeffs,
        bin_width=0.01,
        with_tqdm=True,
        norm=None,
    ) -> tuple[np.ndarray, np.ndarray, float]:
        """
        Broaden dirac delta peaks.

        Parameters
        ----------
            dirac_peaks : np.ndarray
                list of dirac delta peaks
            coeffs : np.ndarray
                Dirac peak intensities
            bin_width : float
                width of the bins to use
            with_tqdm : bool
                whether to use tqdm to show progress of broadening
            norm : float
                value to normalise the spectrum to

        Returns
        -------
            domain : np.ndarray
                range of x-values to used in the spectrum
            data : np.ndarray
                broadened spectrum
            norm_val : float
                value used to normalise the spectrum
        """
        domain = np.arange(self.start, self.stop, bin_width)
        data = np.zeros([len(domain)])
        k_edge_last_x = 0.0

        if coeffs is None:
            coeffs = np.zeros(len(dirac_peaks))

        # Find the peaks in the different broadening regions
        sigma = np.zeros((len(dirac_peaks)))
        mixing = np.zeros((len(dirac_peaks)))

        sigma_precalc = (self.omega_2 - self.omega_1) / (self.ewid_2 - self.ewid_1)
        mix_precalc = (self.mix_2 - self.mix_1) / (self.ewid_2 - self.ewid_1)

        for i in range(len(dirac_peaks)):
            if dirac_peaks[i] <= self.ewid_1:
                sigma[i] = self.omega_1
                mixing[i] = self.mix_1
                # TODO find a better way of doing this
                if dirac_peaks[i] <= self.ewid_1 - 3.0 and dirac_peaks[i] >= 1.0:
                    k_edge_last_x = dirac_peaks[i]
            elif dirac_peaks[i] > self.ewid_2:
                sigma[i] = self.omega_2
                mixing[i] = self.mix_2
            else:
                sigma[i] = self.omega_1 + (
                    sigma_precalc * (dirac_peaks[i] - self.ewid_1)
                )
                mixing[i] = self.mix_1 + (mix_precalc * (dirac_peaks[i] - self.ewid_1))

        if with_tqdm:
            # Parallelise in an openmp style
            if self.n_procs > 1:
                mp_func = partial(
                    Nexafs._mp_broaden,
                    mixing=mixing,
                    dirac_peaks=dirac_peaks,
                    sigma=sigma,
                    coeffs=coeffs,
                )

                # Not able to use tqdm with multiprocessing
                with Pool(self.n_procs) as pool:
                    data = np.array(pool.map(mp_func, domain))

            else:
                data = Nexafs._sp_broaden(
                    data, domain, mixing, dirac_peaks, sigma, coeffs
                )

        else:
            data = Nexafs._sp_broaden(data, domain, mixing, dirac_peaks, sigma, coeffs)

        # Normalise the spectrum
        if norm is None:
            data, norm_val = Nexafs._normalise(data, domain, k_edge_last_x)
            # norm_val = 1
            # print("Not normalising")
        else:
            data, norm_val = Nexafs._normalise(
                data, domain, k_edge_last_x, norm_val=norm
            )
            # norm_val = 1
            # print("Not normalising")

        return domain, data, norm_val

    @staticmethod
    def _plot(angle, path, label, colour=None) -> tuple[np.ndarray, np.ndarray]:
        """
        Private method to parse NEXAFS data and add to mpl graph.

        Parameters
        ----------
            angle : str
                angle of the spectrum
            path : str
                path to the data
            label : str
                label to use for theplot
            colour : str
                colour of the plotted line

        Returns
        -------
            x : np.ndarray
                energy values of broadened spectrum
            y : np.ndarray
                intensity values of broadened spectrum
        """
        x, y = np.loadtxt(
            f"{path}/graphene_Cu_spectrum_{angle}.txt", usecols=(0, 1), unpack=True
        )

        plt.plot(x, y, c=colour, label=label)

        return x, y

    def atom_contribution_plot(
        self,
        deltas: List[str],
        cmap="plasma",
        reversed_cmap=True,
        lower_cmap_range=0.1,
        upper_cmap_range=1.0,
        mpl_figsize: Tuple[float, float] = (8.0, 5.0),
    ) -> None:
        """
        Save a figure of the broadened spectra with individual atom contributions.

        Parameters
        ----------
            deltas : List[str]
                list of x-ray incidence angles
            cmap : str
                colour map to use for invdividual atom contributions in the plot
            reversed_cmap : bool
                whether to reverse the colour map
            lower_cmap_range : float
                normalised lower bound colour of the colour map
            upper_cmap_range : float
                normalised upper bound colour of the colour map
            mpl_figsize : Tuple[float, float]
                size of the saved figure
        """
        # Set the saved figure size
        plt.rcParams["figure.figsize"] = mpl_figsize

        # Reverse the colour map
        if reversed_cmap:
            cmap = plt.get_cmap(cmap).reversed()
        else:
            cmap = plt.get_cmap(cmap)

        # Plot the spectra for all the theta and phi angles
        for i in deltas:
            for j, c in zip(
                self.dirs,
                np.linspace(lower_cmap_range, upper_cmap_range, len(self.dirs)),
            ):
                self._plot(i, f"{j}{i}", j[:-1], cmap(c))

            x, _ = self._plot(i, "./", "Total Spectrum", "black")
            plt.xticks(np.arange(min(x), max(x) + 1, 2))
            plt.ylabel("Normalised Intensity")
            plt.xlabel("Energy (eV)")
            plt.tight_layout()
            plt.legend()

            plt.savefig(f"spectrum_{i}.png")
            print(f"Plotted {i}")

            plt.figure().clear()

    def multi_angle_plot(
        self,
        deltas: List[str],
        cmap="plasma",
        reversed_cmap=False,
        lower_cmap_range=0.1,
        upper_cmap_range=1.0,
        mpl_figsize: Tuple[float, float] = (8.0, 5.0),
    ) -> None:
        """
        Save a figure of the broadened spectra overlayed with spectra angles.

        Parameters
        ----------
            deltas : List[str]
                list of x-ray incidence angles
            cmap : str
                colour map to use for different angles in the plot
            reversed_cmap : bool
                whether to reverse the colour map
            lower_cmap_range : float
                normalised lower bound colour of the colour map
            upper_cmap_range : float
                normalised upper bound colour of the colour map
            mpl_figsize : Tuple[float, float]
                size of the saved figure
        """

        plt.rcParams["figure.figsize"] = mpl_figsize
        x_vals = np.array([])
        x = np.array([])  # Ensure x is bound

        # Reverse the colour map
        if reversed_cmap:
            cmap = plt.get_cmap(cmap).reversed()
        else:
            cmap = plt.get_cmap(cmap)

        for delta, c in zip(
            deltas, np.linspace(lower_cmap_range, upper_cmap_range, len(deltas))
        ):
            x, _ = Nexafs._plot(delta, "./", rf"$\theta = {delta[1:3]}$", cmap(c))
            x_vals = np.concatenate((x_vals, x))

        plt.xticks(np.arange(min(x), max(x) + 1, 2))
        plt.ylabel("Intensity")
        plt.xlabel("Energy (eV)")
        plt.tight_layout()
        plt.legend()

        plt.savefig("spectrum_angles.png")
        print("Plotted angles")


def main(
    n_procs,
    begin,
    end,
    initial_peak,
    omega_1,
    omega_2,
    mix_1,
    mix_2,
    molecule,
    surface,
    excited_atom,
    index_start,
    index_end,
    spectrum_type,
    phi,
    theta,
    excited_nth_element,
    get_i_atoms,
    atom_contribution_plot,
    multi_angle_plot,
):
    nexafs = Nexafs(
        n_procs,
        begin,
        end,
        initial_peak,
        omega_1,
        omega_2,
        mix_1,
        mix_2,
        molecule,
        surface,
        excited_atom,
        index_start,
        index_end,
        spectrum_type,
        phi,
        theta,
        excited_nth_element,
        get_i_atoms,
    )

    # Plot spectrum for all theta and phi angles
    for t in theta:
        for p in phi:
            flat_peaks, flat_bands, peaks, bands = nexafs.get_nexafs_data(
                get_i_atom=get_i_atoms
            )

            print(f"Broadening delta peaks for theta={t} phi={p}...")
            print("Broadening total spectrum...")
            x, y, norm_val = nexafs.broaden(flat_peaks, flat_bands)

            # Write out spectrum to a text file
            np.savetxt(
                f"{molecule}_{surface}_spectrum_t{t}_p{p}.txt", np.vstack((x, y)).T
            )

            if get_i_atoms:
                print("Broadening individual atom spectra...")

                for i in range(len(nexafs.dirs)):
                    print(f"Broadening atom {i+1}...")
                    x, y, _ = nexafs.broaden(
                        peaks[i, :],
                        bands[i, :],
                        with_tqdm=False,
                        norm=norm_val,
                    )

                    # Write out spectrum to a text file
                    np.savetxt(
                        f"{nexafs.dirs[i]}t{t}_p{p}/{molecule}_{surface}_spectrum_t{t}_p{p}.txt",
                        np.vstack((x, y)).T,
                    )

                print("Finished broadening individual atom spectra...")

    # Plot individual atom contributions
    if atom_contribution_plot:
        nexafs.atom_contribution_plot([f"t{t}_p{p}" for t in theta for p in phi])

    # Plot total spectrum for all angles
    if multi_angle_plot:
        nexafs.multi_angle_plot([f"t{t}_p{p}" for t in theta for p in phi])


@click.command()
@click.option(
    "-n",
    "--n_procs",
    default=1,
    type=int,
    show_default=True,
    help="number of processors to use",
)
@click.option(
    "-b", "--begin", default=285.0, type=float, help="start of energy range to plot"
)
@click.option(
    "-e", "--end", default=300.0, type=float, help="end of energy range to plot"
)
@click.option(
    "-i",
    "--initial_peak",
    default=290.0,
    type=float,
    show_default=True,
    help="energy value of the initial peak",
)
@click.option(
    "-o1",
    "--omega_1",
    default=0.75,
    type=float,
    show_default=True,
    help="full width at half maximum for region 1",
)
@click.option(
    "-o2",
    "--omega_2",
    default=2.0,
    type=float,
    show_default=True,
    help="full width at half maximum for region 2",
)
@click.option(
    "-m1",
    "--mix_1",
    default=0.2,
    type=float,
    show_default=True,
    help="Gaussian-Lorentzian mixing parameter for region 1",
)
@click.option(
    "-m2",
    "--mix_2",
    default=0.8,
    type=float,
    show_default=True,
    help="Gaussian-Lorentzian mixing parameter for region 2",
)
@click.option("-m", "--molecule", required=True, type=str, help="adsorbate name")
@click.option("-s", "--surface", required=True, type=str, help="surface element")
@click.option(
    "-a",
    "--excited_atom",
    required=True,
    type=str,
    help="element symbol of the excited atom",
)
@click.option(
    "-i1",
    "--index_start",
    required=True,
    type=int,
    help="index of the first excited atom",
)
@click.option(
    "-i2", "--index_end", required=True, type=int, help="index of the last excited atom"
)
@click.option(
    "-y",
    "--spectrum_type",
    default="avg_polarised",
    type=click.Choice(["tot_summed", "angular", "polarised", "avg_polarised"]),
    show_default=True,
    help="type of spectrum to plot",
)
@click.option(
    "-p", "--phi", default=[60], type=list, show_default=True, help="phi angles"
)
@click.option(
    "-t",
    "--theta",
    default=[0, 25, 53, 90],
    type=list,
    show_default=True,
    help="theta angles",
)
@click.option(
    "-x",
    "--excited_nth_element",
    default=3,
    type=int,
    show_default=True,
    help="the element index of the excited atom of all the elements in the system - this"
    " will always be the last element in the system, so if system contains C, Cu, C:exc"
    " it will be 3",
)
@click.option(
    "-g",
    "--get_i_atoms",
    is_flag=True,
    default=False,
    type=bool,
    show_default=True,
    help="parse and broaden the individual atom spectra",
)
@click.option(
    "-c",
    "--atom_contribution_plot",
    is_flag=True,
    default=True,
    type=bool,
    show_default=True,
    help="plot atom contributions to the total spectrum",
)
@click.option(
    "-l",
    "--multi_angle_plot",
    is_flag=True,
    default=True,
    type=bool,
    show_default=True,
    help="plot the total spectrum for all angles",
)
def nexafs(
    n_procs,
    begin,
    end,
    initial_peak,
    omega_1,
    omega_2,
    mix_1,
    mix_2,
    molecule,
    surface,
    excited_atom,
    index_start,
    index_end,
    spectrum_type,
    phi,
    theta,
    excited_nth_element,
    get_i_atoms,
    atom_contribution_plot,
    multi_angle_plot,
):
    """
    A tool for simulating NEXAFS spectra for surfaces from CASTEP calculations.

    Copyright \u00A9 2023-2024, Dylan Morgan dylan.morgan@warwick.ac.uk
    """
    main(
        n_procs,
        begin,
        end,
        initial_peak,
        omega_1,
        omega_2,
        mix_1,
        mix_2,
        molecule,
        surface,
        excited_atom,
        index_start,
        index_end,
        spectrum_type,
        phi,
        theta,
        excited_nth_element,
        get_i_atoms,
        atom_contribution_plot,
        multi_angle_plot,
    )


if __name__ == "__main__":
    nexafs()
