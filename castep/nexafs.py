#!/usr/bin/env python3

import os
from dataclasses import dataclass
from functools import partial
from multiprocessing import Pool
from typing import Annotated, Literal

import click
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from tqdm import tqdm


@dataclass
class GreaterThan:
    """
    Represent a value greater than a given value.
    """

    min: float


class NEXAFS:
    def __init__(
        self,
        n_procs: int,
        begin: float,
        end: float,
        initial_peak: float,
        omega_1: float,
        omega_2: float,
        mix_1: float,
        mix_2: float,
        molecule: str,
        surface: str,
        excited_atom: str,
        index_start: int,
        index_end: int,
        spectrum_type: Literal["tot_summed", "angular", "polarised", "avg_polarised"],
        phi: list[str],
        theta: list[str],
        excited_nth_element: int,
        plot_i_atoms: bool,
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
        spectrum_type : Literal["tot_summed", "angular", "polarised", "avg_polarised"]
            type of spectrum to plot
        phi : list[str]
            phi angles
        theta : list[str]
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
        spectrum_type_opts = ["tot_summed", "angular", "polarised", "avg_polarised"]
        self.n_type = spectrum_type_opts.index(spectrum_type) + 1

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
            f"{self.element}{str(atom_ids[0])}/t{theta[0]}_p{phi[0]}{self.fname}"
        )
        tmp_bands_l = len(tmp_bands)

        # Create arrays with sizes of the system to use
        self.peaks = np.zeros([len(atom_ids), tmp_bands_l])
        self.bands = np.zeros([len(atom_ids), tmp_bands_l])

    @property
    def k_edge_max(self) -> np.float64:
        """
        Get the maximum intensity of the first peak.

        Returns
        -------
        np.float64
            maximum intensity of the first peak
        """

        return np.max(self.bands[np.where(self.peaks <= self.ewid_1)])

    def check_prev_broadening(self, theta: str, phi: str) -> bool:
        """
        Check if broadened deltas files already exist.

        Returns
        -------
            True if broadened deltas files exist, False otherwise
        """

        if os.path.isfile(f"{self.molecule}_{self.metal}_deltas_t{theta}_p{phi}.txt"):
            return True
        else:
            return False

    def get_exp_nexafs_data(self, path: str) -> tuple[np.ndarray, np.ndarray]:
        """
        Parse the experimental NEXAFS data

        Data must be in the format of two columns: energy and intensity, not peaks and
        bands.

        Parameters
        ----------
        path : str
            path to the NEXAFS data file

        Returns
        -------
        x : npt.NDArray[np.float64]
            energy values of the broadened spectrum
        y : npt.NDArray[np.float64]
            intensity values of the broadened spectrum
        """

        x, y = np.loadtxt(path, unpack=True)

        return x, y

    def get_sim_nexafs_data(
        self,
        theta: str,
        phi: str,
        get_i_atom: bool = False,
    ) -> tuple[
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
    ]:
        """
        Parse the simulated NEXAFS data for all given theta and phi angles.

        Parameters
        ----------
        peaks : npt.NDArray[np.float64]
            array to store the Delta peak positions
        bands : npt.NDArray[np.float64]
            array to store the Delta peak intensities
        get_i_atom : bool, default=False
            whether to parse and save the individual atom spectra to a file

        Returns
        -------
        flat_peaks : npt.NDArray[np.float64]
            flattened array of Delta peak positions
        flat_bands : npt.NDArray[np.float64]
            flattened array of Delta peak intensities
        peaks : npt.NDArray[np.float64]
            array of Delta peak positions
        bands : npt.NDArray[np.float64]
            array of Delta peak intensities
        """

        print(f"Parsing NEXAFS data for theta={theta} phi={phi}...")
        for i in range(len(self.dirs)):
            nexafs = np.loadtxt(
                f"{self.dirs[i]}t{theta}_p{phi}{self.fname}", skiprows=1
            )
            self.peaks[i, :] = nexafs[:, 0]
            self.bands[i, :] = nexafs[:, self.n_type]

            if get_i_atom:  # Save each atom spectrum to a file
                i_atom_peaks = nexafs[:, 0]
                i_atom_bands = nexafs[:, self.n_type]
                np.savetxt(
                    f"{self.dirs[i]}t{theta}_p{phi}/{self.molecule}_"
                    f"{self.metal}_deltas_t{theta}_p{phi}.txt",
                    np.vstack((i_atom_peaks, i_atom_bands)).T,
                    header="# <x in eV> Intensity",
                )

        flat_peaks = self.peaks.flatten()
        flat_bands = self.bands.flatten()

        # Write out all of the data into a delta peaks file
        np.savetxt(
            f"{self.molecule}_{self.metal}_deltas_t{theta}_p{phi}.txt",
            np.vstack((flat_peaks, flat_bands)).T,
            header="# <x in eV> Intensity",
        )

        print(f"Finished writing out delta peaks file for theta={theta} phi={phi}")

        return flat_peaks, flat_bands, self.peaks, self.bands

    def rigid_shift(self, shift_val: float):
        """
        Rigidly shift the spectrum by a value in the x-direction

        Parameters
        ----------
        shift_val : float
            amount to shift the spectrum by
        """

    # def normalise(self, data, domain, ewid_1) -> tuple[npt.NDArray[np.float64], float]:
    #     """
    #     Normalise spectrum such that first peak has intensity of 1.

    #     Parameters
    #     ----------
    #     data : np.ndarray
    #         broadened spectrum to normalise
    #     domain : np.ndarray
    #         range of x-values to include in the spectrum
    #     ewid_1 : float
    #         x-value at which the first peak region ends

    #     Returns
    #     -------
    #     data : np.ndarray
    #         normalised spectrum
    #     k_edge_max : float
    #         maximum intensity of the first peak
    #     """

    #     k_edge_max = np.max(data[np.where(domain <= ewid_1)])

    #     # Scale the rest of the data proportionally
    #     data /= k_edge_max

    #     return data, k_edge_max

    # TODO: implement @lru_cache(maxsize=128)
    @staticmethod
    def _schmid_pseudo_voigt(
        domain: npt.NDArray[np.float64],
        m: float,
        E: float,
        omega: Annotated[float, GreaterThan(0)],
    ) -> npt.NDArray[np.float64]:
        """
        Apply broadening scheme for XPS spectra. See the following reference for more details:
        https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/sia.5521

        Parameters
        ----------
        domain : npt.NDArray[np.float64]
            range of x-values to include in the spectrum
        m : float
            Gaussian-Lorentzian mixing parameter
        E : float
            line centre/delta peak
        omega : Annotated[float, GreaterThan(0)]
            full width at half maximum (omega > 0)

        Returns
        -------
        data : npt.NDArray[np.float64]
            broadened spectrum
        """

        return (1 - m) * np.sqrt((4 * np.log(2)) / (np.pi * omega**2)) * np.exp(
            -(4 * np.log(2) / omega**2) * (domain - E) ** 2
        ) + m * (1 / (2 * np.pi)) * (omega / ((omega / 2) ** 2 + (domain - E) ** 2))

    @staticmethod
    def _mp_broaden(
        domain, mixing, delta_peaks, sigma, coeffs
    ) -> npt.NDArray[np.float64]:
        """
        Perform the broadening in parallel.

        Parameters
        ----------
        domain : np.ndarray
            range of x-values to include in the spectrum
        mixing : np.ndarray
            Gaussian-Lorentzian mixing parameter
        delta_peaks : np.ndarray
            line centre/delta peak
        sigma : np.ndarray
            linearly interpolated full width at half maximum
        coeffs : np.ndarray
            Delta peak intensities

        Returns
        -------
        data : np.ndarray
            broadened spectrum
        """

        data = np.sum(
            NEXAFS._schmid_pseudo_voigt(domain, mixing, delta_peaks, sigma) * coeffs
        )

        return data

    @staticmethod
    def _sp_broaden(
        data, domain, mixing, delta_peaks, sigma, coeffs
    ) -> npt.NDArray[np.float64]:
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
        delta_peaks : np.ndarray
            line centre/delta peak
        sigma : np.ndarray
            linearly interpolated full width at half maximum
        coeffs : np.ndarray
            Delta peak intensities

        Returns
        -------
        data : np.ndarray
            broadened spectrum
        """

        for i in tqdm(range(len(domain))):
            data[i] = np.sum(
                NEXAFS._schmid_pseudo_voigt(domain[i], mixing, delta_peaks, sigma)
                * coeffs
            )

        return data

    def broaden(
        self,
        delta_peaks: npt.NDArray[np.float64],
        coeffs: npt.NDArray[np.float64],
        bin_width=0.01,
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """
        Broaden delta delta peaks.

        Parameters
        ----------
        delta_peaks : npt.NDArray[np.float64]
            list of delta delta peaks
        coeffs : npt.NDArray[np.float64]
            Delta peak intensities
        bin_width : float, default=0.01
            width of the bins to use

        Returns
        -------
        domain : npt.NDArray[np.float64]
            range of x-values to used in the spectrum
        data : npt.NDArray[np.float64]
            broadened spectrum
        """

        domain = np.arange(self.start, self.stop, bin_width)
        data = np.zeros([len(domain)])
        # k_edge_last_x = 0.0

        if coeffs is None:
            coeffs = np.zeros(len(delta_peaks))

        # Find the peaks in the different broadening regions
        sigma = np.zeros((len(delta_peaks)))
        mixing = np.zeros((len(delta_peaks)))

        sigma_precalc = (self.omega_2 - self.omega_1) / (self.ewid_2 - self.ewid_1)
        mix_precalc = (self.mix_2 - self.mix_1) / (self.ewid_2 - self.ewid_1)

        for i in range(len(delta_peaks)):
            if delta_peaks[i] <= self.ewid_1:
                sigma[i] = self.omega_1
                mixing[i] = self.mix_1
                # TODO find a better way of doing this
                # if delta_peaks[i] <= self.ewid_1 - 3.0 and delta_peaks[i] >= 1.0:
                #     k_edge_last_x = delta_peaks[i]
            elif delta_peaks[i] > self.ewid_2:
                sigma[i] = self.omega_2
                mixing[i] = self.mix_2
            else:
                sigma[i] = self.omega_1 + (
                    sigma_precalc * (delta_peaks[i] - self.ewid_1)
                )
                mixing[i] = self.mix_1 + (mix_precalc * (delta_peaks[i] - self.ewid_1))

        # Parallelise in an openmp style
        if self.n_procs > 1:
            mp_func = partial(
                NEXAFS._mp_broaden,
                mixing=mixing,
                delta_peaks=delta_peaks,
                sigma=sigma,
                coeffs=coeffs,
            )

            with Pool(self.n_procs) as pool:
                data = np.array(pool.map(mp_func, domain))

        else:
            data = NEXAFS._sp_broaden(data, domain, mixing, delta_peaks, sigma, coeffs)

        return domain, data

    @staticmethod
    def _plot(
        angle, path, label, colour=None
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
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
        colour : str, default=None
            colour of the plotted line

        Returns
        -------
        x : npt.NDArray[np.float64]
            energy values of broadened spectrum
        y : npt.NDArray[np.float64]
            intensity values of broadened spectrum
        """

        x, y = np.loadtxt(
            f"{path}/graphene_Cu_spectrum_{angle}.txt", usecols=(0, 1), unpack=True
        )

        plt.plot(x, y, c=colour, label=label)

        return x, y

    def atom_contribution_plot(
        self,
        deltas: list[str],
        cmap="plasma",
        reversed_cmap=True,
        lower_cmap_range=0.1,
        upper_cmap_range=1.0,
        mpl_figsize: tuple[float, float] = (8.0, 5.0),
    ) -> None:
        """
        Save a figure of the broadened spectra with individual atom contributions.

        Parameters
        ----------
        deltas : list[str]
            list of x-ray incidence angles
        cmap : str, default="plasma"
            colour map to use for invdividual atom contributions in the plot
        reversed_cmap : bool, default=True
            whether to reverse the colour map
        lower_cmap_range : float, default=0.1
            normalised lower bound colour of the colour map
        upper_cmap_range : float, default=1.0
            normalised upper bound colour of the colour map
        mpl_figsize : Tuple[float, float], default=(8.0, 5.0)
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
            # plt.ylabel("Normalised Intensity")
            plt.ylabel("Intensity")
            plt.xlabel("Energy (eV)")
            plt.tight_layout()
            plt.legend()

            plt.savefig(f"spectrum_{i}.png")
            print(f"Plotted {i}")

            plt.figure().clear()

    def multi_angle_plot(
        self,
        deltas: list[str],
        cmap="plasma",
        reversed_cmap=False,
        lower_cmap_range=0.1,
        upper_cmap_range=1.0,
        mpl_figsize: tuple[float, float] = (8.0, 5.0),
    ) -> None:
        """
        Save a figure of the broadened spectra overlayed with spectra angles.

        Parameters
        ----------
            deltas : list[str]
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
            x, _ = NEXAFS._plot(delta, "./", rf"$\theta = {delta[1:3]}$", cmap(c))
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
    nexafs = NEXAFS(
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
            # Check if deltas files already exist
            # if nexafs.check_prev_broadening(t, p):
            #     print(f"Delta peaks file for theta={t} phi={p} found, skipping...")
            #     continue

            flat_peaks, flat_bands, peaks, bands = nexafs.get_sim_nexafs_data(
                t, p, get_i_atom=get_i_atoms
            )

            print(f"Broadening delta peaks for theta={t} phi={p}...")
            print("Broadening total spectrum...")
            x, y = nexafs.broaden(flat_peaks, flat_bands)

            # Write out spectrum to a text file
            np.savetxt(
                f"{molecule}_{surface}_spectrum_t{t}_p{p}.txt", np.vstack((x, y)).T
            )

            if get_i_atoms:
                print("Broadening individual atom spectra...")

                for i in tqdm(range(len(nexafs.dirs))):
                    x, y = nexafs.broaden(peaks[i, :], bands[i, :])

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
    type=click.IntRange(1, clamp=True),
    show_default=True,
    help="number of processors to use",
)
@click.option(
    "-b", "--begin", required=True, type=float, help="start of energy range to plot"
)
@click.option(
    "-e", "--end", required=True, type=float, help="end of energy range to plot"
)
@click.option(
    "-i",
    "--initial_peak",
    required=True,
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
    "-p", "--phi", default=["60"], type=list[str], show_default=True, help="phi angles"
)
@click.option(
    "-t",
    "--theta",
    default=["00", "25", "53", "90"],
    type=list[str],
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
    show_default=True,
    help="parse and broaden the individual atom spectra",
)
@click.option(
    "-c",
    "--atom_contribution_plot",
    is_flag=True,
    default=False,
    show_default=True,
    help="plot atom contributions to the total spectrum",
)
@click.option(
    "-l",
    "--multi_angle_plot",
    is_flag=True,
    default=False,
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
