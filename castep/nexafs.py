#!/usr/bin/env python3

import os
from copy import copy
from dataclasses import dataclass
from functools import partial
from multiprocessing import Pool
from typing import Annotated, Literal

import click
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
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

    # Set the x axis of the spectrum to be in integer 2 eV increments
    plt.figure().gca().xaxis.set_major_locator(ticker.MultipleLocator(2))

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
        root_dir: str = "./",
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
        root_dir : str, default="./"
            root directory of the data
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

        # Set up a list of all the directories all the data is in
        self.root_dir = root_dir
        self.dirs = np.array([])
        for n in atom_ids:
            direc = f"{self.root_dir}/{excited_atom}{str(n)}/"
            if os.path.isdir(direc):
                self.dirs = np.append(self.dirs, direc)

        # Don't attempt to parse if not foound ie. for experimental data
        if len(self.dirs) > 0:
            # Deltas file for each molpdos run
            self.fname = f"{self.molecule}_{self.metal}_{self.atom}_1_1_1_deltas.dat"

            # Get the length of the deltas file
            num_peaks = len(
                np.loadtxt(f"{self.dirs[0]}/t{theta[0]}_p{phi[0]}/{self.fname}")
            )

            # Create arrays with sizes of the system to use
            self.separate_peaks = np.zeros([len(self.dirs), num_peaks])
            self.separate_bands = np.zeros([len(self.dirs), num_peaks])

        else:
            self.separate_peaks = None
            self.separate_bands = None

        self.total_peaks = None
        self.total_bands = None

        self.broadened = False

    @property
    def separate_peaks(self):
        """Peaks where each row is an atom contributing to the spectrum."""

        return self._separate_peaks

    @separate_peaks.setter
    def separate_peaks(self, peaks):
        self._separate_peaks = peaks

    @property
    def total_peaks(self):
        """Flattened peaks array."""

        return self._total_peaks

    @total_peaks.setter
    def total_peaks(self, peaks):
        self._total_peaks = peaks

    @property
    def separate_bands(self):
        """Bands where each row is an atom contributing to the spectrum."""

        return self._separate_bands

    @separate_bands.setter
    def separate_bands(self, bands):
        self._separate_bands = bands

    @property
    def total_bands(self):
        """Flattened bands array."""

        return self._total_bands

    @total_bands.setter
    def total_bands(self, bands):
        self._total_bands = bands

    @property
    def broadened(self):
        """Whether the spectrum has been broadened."""

        return self._broadened

    @broadened.setter
    def broadened(self, broadened: bool):
        self._broadened = broadened

    @property
    def k_edge_max(self):
        """The maximum intensity of the broadened first peak."""

        if self.total_bands is None:
            self._k_edge_max = np.max(
                self.separate_bands[
                    np.asarray(self.separate_peaks < self.ewid_1).nonzero()
                ]
            )

        else:
            self._k_edge_max = np.max(
                self.total_bands[np.asarray(self.total_peaks < self.ewid_1).nonzero()]
            )

        return self._k_edge_max

    @property
    def energy_k_edge_max(self):
        """The energy of the maximum intensity of the broadened first peak."""

        if not self.broadened:
            # Don't allow calculation of this if the broadened peaks have not been calculated
            raise ValueError(
                f"{NEXAFS.energy_k_edge_max.fget.__name__}() can only be called after the spectrum has been broadened."  # pyright: ignore
            )

        else:
            return self.total_peaks[
                np.argmax(np.asarray(self.total_peaks < self.ewid_1).nonzero())
            ]

    def parse_experimental(self, file_name: str) -> None:
        """
        Parse the experimental NEXAFS data

        Data must be in the format of two columns: energy and intensity, not peaks and
        bands.

        Parameters
        ----------
        file_name : str
            name of the file to parse containing the experimental data
        """

        self.broadened = True

        self.total_peaks, self.total_bands = np.loadtxt(
            f"{self.root_dir}/{file_name}", unpack=True
        )

    def parse_simulated(self, theta: str, phi: str, get_i_atom: bool = False) -> None:
        """
        Parse the simulated NEXAFS data for all given theta and phi angles.

        Parameters
        ----------
        theta : str
            theta angle of incidence X-ray to parse
        phi : str
            phi angle of incidence X-ray to parse
        get_i_atom : bool, default=False
            whether to parse individual atom spectra
        """

        print(f"Parsing NEXAFS data for theta={theta} phi={phi}...")
        for i in range(len(self.dirs)):
            nexafs = np.loadtxt(
                f"{self.dirs[i]}t{theta}_p{phi}/{self.fname}",
                skiprows=1,
                unpack=True,
            )
            self.separate_peaks[i] = nexafs[0]
            self.separate_bands[i] = nexafs[self.n_type]

            if get_i_atom:  # Save each atom spectrum to a file
                i_atom_peaks = nexafs[0]
                i_atom_bands = nexafs[self.n_type]
                np.savetxt(
                    f"{self.dirs[i]}t{theta}_p{phi}/{self.molecule}_"
                    f"{self.metal}_deltas_t{theta}_p{phi}.txt",
                    np.vstack((i_atom_peaks, i_atom_bands)).T,
                    header="# <x in eV> Intensity",
                )

        self.total_peaks = self.separate_peaks.flatten()
        self.total_bands = self.separate_bands.flatten()

        # Write out all of the data into a delta peaks file
        np.savetxt(
            f"{self.root_dir}/{self.molecule}_{self.metal}_deltas_t{theta}_p{phi}.txt",
            np.vstack((self.total_peaks, self.total_bands)).T,
            header="# <x in eV> Intensity",
        )

        print(f"Finished writing out delta peaks file for theta={theta} phi={phi}")

    @staticmethod
    def _schmid_pseudo_voigt(
        domain: npt.NDArray[np.float64],
        m: npt.NDArray[np.float64],
        E: float,
        omega: npt.NDArray[Annotated[np.float64, GreaterThan(0)]],
    ) -> npt.NDArray[np.float64]:
        """
        Apply broadening scheme for XPS spectra. See the following reference for more details:
        https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/sia.5521

        Parameters
        ----------
        domain : npt.NDArray[np.float64]
            range of x-values to include in the spectrum
        m : npt.NDArray[np.float64]
            Gaussian-Lorentzian mixing parameter
        E : float
            line centre/delta peak
        omega : npt.NDArray[Annnotated[np.float64, GreaterThan(0)]
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
        domain: npt.NDArray[np.float64],
        mixing: npt.NDArray[np.float64],
        delta_peak: float,
        sigma: npt.NDArray[np.float64],
        coeffs: npt.NDArray[np.float64],
    ) -> np.float64:
        """
        Perform the broadening in parallel.

        Parameters
        ----------
        domain : npt.NDArray[np.float64]
            range of x-values to include in the spectrum
        mixing : npt.NDArray[np.float64]
            Gaussian-Lorentzian mixing parameter
        delta_peak : float
            line centre/delta peak
        sigma : npt.NDArray[np.float64]
            linearly interpolated full width at half maximum
        coeffs : npt.NDArray[np.float64]
            Delta peak intensities

        Returns
        -------
        data : float64
            broadened spectrum
        """

        data = np.sum(
            NEXAFS._schmid_pseudo_voigt(domain, mixing, delta_peak, sigma) * coeffs
        )

        return data

    @staticmethod
    def _sp_broaden(
        data: npt.NDArray[np.float64],
        domain: npt.NDArray[np.float64],
        mixing: npt.NDArray[np.float64],
        delta_peak: float,
        sigma: npt.NDArray[np.float64],
        coeffs: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """
        Perform the broadening in serial.

        Parameters
        ----------
        data : npt.NDArray[np.float64]
            non-broadened spectrum
        domain : npt.NDArray[np.float64]
            range of x-values to include in the spectrum
        mixing : npt.NDArray[np.float64]
            Gaussian-Lorentzian mixing parameter
        delta_peak : float
            line centre/delta peak
        sigma : npt.NDArray[np.float64]
            linearly interpolated full width at half maximum
        coeffs : npt.NDArray[np.float64]
            Delta peak intensities

        Returns
        -------
        data : npt.NDArray[np.float64]
            broadened spectrum
        """

        for i in tqdm(range(len(domain))):
            data[i] = np.sum(
                NEXAFS._schmid_pseudo_voigt(domain[i], mixing, delta_peak, sigma)
                * coeffs
            )

        return data

    def broaden(self, bin_width=0.01) -> None:
        """
        Broaden NEXAFS delta peaks.

        Parameters
        ----------
        bin_width : float, default=0.01
            width of the bins to use

        Returns
        -------
        domain : npt.NDArray[np.float64]
            range of x-values to use in the spectrum
        data : npt.NDArray[np.float64]
            broadened spectrum
        """

        domain = np.arange(self.start, self.stop, bin_width)
        broadened = np.zeros([len(domain)])

        if self.total_bands is None:
            coeffs = np.zeros(len(self.total_peaks))
        else:
            coeffs = self.total_bands

        # Find the peaks in the different broadening regions
        sigma = np.zeros((len(self.total_peaks)))
        mixing = np.zeros((len(self.total_peaks)))

        sigma_precalc = (self.omega_2 - self.omega_1) / (self.ewid_2 - self.ewid_1)
        mix_precalc = (self.mix_2 - self.mix_1) / (self.ewid_2 - self.ewid_1)

        for i in range(len(self.total_peaks)):
            if self.total_peaks[i] <= self.ewid_1:
                sigma[i] = self.omega_1
                mixing[i] = self.mix_1
                # TODO find a better way of doing this
                # if delta_peaks[i] <= self.ewid_1 - 3.0 and delta_peaks[i] >= 1.0:
                #     k_edge_last_x = delta_peaks[i]
            elif self.total_peaks[i] > self.ewid_2:
                sigma[i] = self.omega_2
                mixing[i] = self.mix_2
            else:
                sigma[i] = self.omega_1 + (
                    sigma_precalc * (self.total_peaks[i] - self.ewid_1)
                )
                mixing[i] = self.mix_1 + (
                    mix_precalc * (self.total_peaks[i] - self.ewid_1)
                )

        # Parallelise in an openmp style
        if self.n_procs > 1:
            mp_func = partial(
                NEXAFS._mp_broaden,
                mixing=mixing,
                delta_peak=self.total_peaks,
                sigma=sigma,
                coeffs=coeffs,
            )

            with Pool(self.n_procs) as pool:
                broadened = np.array(pool.map(mp_func, domain))

        else:
            broadened = NEXAFS._sp_broaden(
                broadened, domain, mixing, self.total_peaks, sigma, coeffs
            )

        self.total_peaks = domain
        self.total_bands = broadened
        self.broadened = True

    def rigid_shift(self, shift_val: float) -> None:
        """
        Rigidly shift the spectrum by a value in the x-direction

        Parameters
        ----------
        shift_val : float
            amount to shift the spectrum by
        """

        self.ewid_1 -= shift_val
        self.ewid_2 -= shift_val

        self.start -= shift_val
        self.stop -= shift_val

        self.separate_peaks -= shift_val
        self.total_peaks -= shift_val

    def normalise(self) -> None:
        """Normalise the broadened spectrum so the k-edge has an intensity of 1"""

        if self.separate_bands is not None:
            self.separate_bands /= self.k_edge_max

        if self.total_bands is not None:
            self.total_bands /= self.k_edge_max

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
            label to use for the plot
        colour : str, default=None
            colour of the plotted line

        Returns
        -------
        x : npt.NDArray[np.float64]
            energy values of broadened spectrum
        y : npt.NDArray[np.float64]
            intensity values of broadened spectrum
        """

        x, y = np.loadtxt(f"{path}/graphene_Cu_spectrum_{angle}.txt", unpack=True)

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
    sim_initial_peak,
    exp_initial_peak,
    root_dir,
    experimental_data,
    edge_atom_root,
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

    exp = NEXAFS(
        n_procs,
        begin,
        end,
        exp_initial_peak,
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
        root_dir="/".join(experimental_data.split("/")[:-1]),
    )

    sw_sim = NEXAFS(
        n_procs,
        begin,
        end,
        sim_initial_peak,
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
        root_dir=root_dir,
    )

    edge_atom = NEXAFS(
        n_procs,
        begin,
        end,
        sim_initial_peak,
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
        root_dir=edge_atom_root,
    )

    # Get the experimental data and normalise
    exp.parse_experimental(experimental_data.split("/")[-1])
    exp.normalise()

    # Plot spectrum for all theta and phi angles
    for t in theta:
        for p in phi:
            sw_sim.parse_simulated(t, p, get_i_atom=get_i_atoms)
            edge_atom.parse_simulated(t, p, get_i_atom=get_i_atoms)

            mix = copy(sw_sim)
            mix.separate_bands = np.append(
                mix.separate_bands, edge_atom.total_bands * (16 * 5)
            )
            mix.total_bands = mix.separate_bands.flatten()
            mix.separate_peaks = np.append(mix.separate_peaks, edge_atom.total_peaks)
            mix.total_peaks = mix.separate_peaks.flatten()

            # print(f"Broadening delta peaks for theta={t} phi={p}...")
            print(f"Broadening delta peaks for theta={theta[0]} phi={phi[0]}...")
            print("Broadening total spectrum...")
            edge_atom.broaden()
            edge_atom.normalise()
            edge_atom.rigid_shift(edge_atom.energy_k_edge_max - exp.energy_k_edge_max)
            mix.broaden()
            mix.normalise()
            mix.rigid_shift(mix.energy_k_edge_max - exp.energy_k_edge_max)
            sw_sim.broaden()
            sw_sim.normalise()
            sw_sim.rigid_shift(sw_sim.energy_k_edge_max - exp.energy_k_edge_max)

    if get_i_atoms:
        raise NotImplementedError
        # print("Broadening individual atom spectra...")

        # for i in tqdm(range(len(sw_sim.dirs))):
        #     sw_sim.broaden()

        #     # Write out spectrum to a text file
        #     np.savetxt(
        #         f"{sw_sim.dirs[i]}t{t}_p{p}/{molecule}_{surface}_spectrum_t{t}_p{p}.txt",
        #         np.vstack((sw_sim.total_peaks, sw_sim.total_bands)).T,
        #     )

        # print("Finished broadening individual atom spectra...")

    plt.plot(sw_sim.total_peaks, sw_sim.total_bands, "--", label="Defect Atoms")
    plt.plot(edge_atom.total_peaks, edge_atom.total_bands, "--", label="Edge Atom")
    plt.plot(exp.total_peaks, exp.total_bands, label="Experiment")
    plt.plot(mix.total_peaks, mix.total_bands, label="Mixed")
    plt.xlim(min(sw_sim.total_peaks), max(sw_sim.total_peaks))
    plt.xlabel("Energy / eV")
    plt.ylabel("Normalised Intensity")
    plt.legend()
    plt.savefig("broadened.png")

    # Plot individual atom contributions
    if atom_contribution_plot:
        raise NotImplementedError
        # nexafs.atom_contribution_plot([f"t{t}_p{p}" for t in theta for p in phi])

    # Plot total spectrum for all angles
    if multi_angle_plot:
        raise NotImplementedError
        # nexafs.multi_angle_plot([f"t{t}_p{p}" for t in theta for p in phi])


@click.command()
@click.option(
    "-n",
    "--n_procs",
    default=1,
    type=click.IntRange(1, clamp=True),
    show_default=True,
    help="number of processors to use",
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
    "-b",
    "--begin",
    required=True,
    type=float,
    help="start of energy range to plot for the non-shifted spectrum",
)
@click.option(
    "-e",
    "--end",
    required=True,
    type=float,
    help="end of energy range to plot for the non-shifted spectrum",
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
    "-si",
    "--simulated_initial_peak",
    required=True,
    type=float,
    help="energy value of the simulated initial peak",
)
@click.option(
    "-ei",
    "--experimental_initial_peak",
    type=float,
    help="energy value of the experimental initial peak",
)
@click.option(
    "-r",
    "--root_dir",
    default="./",
    type=click.Path(exists=True, file_okay=False),
    show_default=True,
    help="root directory of the data",
)
@click.option(
    "-ed",
    "--experimental_data",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="path to the experimental data to plot",
)
@click.option(
    "-ar",
    "--edge_atom_root",
    required=True,
    type=click.Path(exists=True, file_okay=False),
    help="path to the edge atom root directory",
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
@click.option(
    "-y",
    "--spectrum_type",
    default="avg_polarised",
    type=click.Choice(["tot_summed", "angular", "polarised", "avg_polarised"]),
    show_default=True,
    help="type of spectrum to plot",
)
@click.option(
    "-p",
    "--phi",
    default=["60"],
    show_default=True,
    multiple=True,
    help="phi angles",
)
@click.option(
    "-t",
    "--theta",
    default=["00", "25", "53", "90"],
    show_default=True,
    multiple=True,
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
    molecule,
    surface,
    excited_atom,
    begin,
    end,
    index_start,
    index_end,
    simulated_initial_peak,
    experimental_initial_peak,
    root_dir,
    experimental_data,
    edge_atom_root,
    omega_1,
    omega_2,
    mix_1,
    mix_2,
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
        simulated_initial_peak,
        experimental_initial_peak,
        root_dir,
        experimental_data,
        edge_atom_root,
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
