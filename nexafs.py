#!/usr/bin/env python3

import os
from functools import partial
from multiprocessing import Pool

import click
import numpy as np
from tqdm import tqdm


class Nexafs:
    def __init__(self):
        pass

    def get_nexafs_data(
        theta, phi, dirs, fname, peaks, bands, n_type, molecule, metal, get_i_atom=False
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Parse the NEXAFS data for all theta and phi"""

        print(f"Parsing NEXAFS data for theta={theta} phi={phi}...")
        for i in range(len(dirs)):
            nexafs = np.loadtxt(f"{dirs[i]}t{theta}_p{phi}{fname}")
            peaks[i, :] = nexafs[:, 0]
            bands[i, :] = nexafs[:, n_type]

            if get_i_atom:  # Save each atom spectrum to a file
                i_atom_peaks = nexafs[:, 0]
                i_atom_bands = nexafs[:, n_type]
                np.savetxt(
                    f"{dirs[i]}t{theta}_p{phi}/{molecule}_{metal}_deltas_t{theta}_p{phi}.txt",
                    np.vstack((i_atom_peaks, i_atom_bands)).T,
                    header="# <x in eV> Intensity",
                )

        flat_peaks = peaks.flatten()
        flat_bands = bands.flatten()

        # Write out all of the data into a delta peaks file
        np.savetxt(
            f"{molecule}_{metal}_deltas_t{theta}_p{phi}.txt",
            np.vstack((flat_peaks, flat_bands)).T,
            header="# <x in eV> Intensity",
        )

        print(f"Finished writing out delta peaks file for theta={theta} phi={phi}")

        return flat_peaks, flat_bands, peaks, bands

    def _schmid_pseudo_voigt(domain, m, E, omega) -> np.ndarray:
        """
        Apply broadening scheme for XPS spectra
        https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/sia.5521

        domain = linspace of x range and bin width
        m = Gaussian-Lorentzian mixing parameter
        E = line centre (aka dirac peak)
        omega = full width at half maximum (omega > 0)
        """

        return (1 - m) * np.sqrt((4 * np.log(2)) / (np.pi * omega**2)) * np.exp(
            -(4 * np.log(2) / omega**2) * (domain - E) ** 2
        ) + m * (1 / (2 * np.pi)) * (omega / ((omega / 2) ** 2 + (domain - E) ** 2))

    def _normalise(data, domain, ewid_1, norm_val=None):
        """Normalise spectrum such that first peak has intensity of 1."""

        # Find height of first peak
        if norm_val is None:
            k_edge_max = np.max(data[np.where(domain <= ewid_1)])
        else:
            k_edge_max = norm_val

        # Scale the rest of the data proportionally
        data /= k_edge_max

        return data, k_edge_max

    def _mp_broaden(domain, mixing, dirac_peaks, sigma, coeffs):
        data = np.sum(_schmid_pseudo_voigt(domain, mixing, dirac_peaks, sigma) * coeffs)

        return data

    def broaden(
        n_procs,
        start,
        stop,
        dirac_peaks,
        coeffs,
        omega_1,
        omega_2,
        mix_1,
        mix_2,
        ewid_1,
        ewid_2,
        bin_width=0.01,
        with_tqdm=True,
        norm=None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Broaden dirac delta peaks

        start = beginning of x range
        stop = end of x range
        dirac_peaks = list of dirac delta peaks
        coeffs =
        omega = full width at half maximum (omega > 0)
        mix = Gaussian-Lorentzian mixing parameter
        ewid = x-values at which to change mix
        eta = full width at half maximum (eta > 0)
        """
        domain = np.arange(start, stop, bin_width)
        data = np.zeros([len(domain)])
        k_edge_last_x = 0.0

        if coeffs is None:
            coeffs = np.zeros(len(dirac_peaks))

        # Find the peaks in the different broadening regions
        sigma = np.zeros((len(dirac_peaks)))
        mixing = np.zeros((len(dirac_peaks)))

        sigma_precalc = (omega_2 - omega_1) / (ewid_2 - ewid_1)
        mix_precalc = (mix_2 - mix_1) / (ewid_2 - ewid_1)

        for i in range(len(dirac_peaks)):
            if dirac_peaks[i] <= ewid_1:
                sigma[i] = omega_1
                mixing[i] = mix_1
                # TODO find a better way of doing this
                if dirac_peaks[i] <= ewid_1 - 3.0 and dirac_peaks[i] >= 1.0:
                    k_edge_last_x = dirac_peaks[i]
            elif dirac_peaks[i] > ewid_2:
                sigma[i] = omega_2
                mixing[i] = mix_2
            else:
                sigma[i] = omega_1 + (sigma_precalc * (dirac_peaks[i] - ewid_1))
                mixing[i] = mix_1 + (mix_precalc * (dirac_peaks[i] - ewid_1))

        if with_tqdm:
            # Parallelise in an openmp style
            if n_procs > 1:
                mp_func = partial(
                    _mp_broaden,
                    mixing=mixing,
                    dirac_peaks=dirac_peaks,
                    sigma=sigma,
                    coeffs=coeffs,
                )

                # Not able to use tqdm with multiprocessing
                with Pool(n_procs) as pool:
                    data = np.array(pool.map(mp_func, domain))

            else:
                for i in tqdm(range(len(domain))):
                    data[i] = np.sum(
                        _schmid_pseudo_voigt(domain[i], mixing, dirac_peaks, sigma)
                        * coeffs
                    )

        else:
            for i in tqdm(range(len(domain))):
                data[i] = np.sum(
                    _schmid_pseudo_voigt(domain[i], mixing, dirac_peaks, sigma) * coeffs
                )

        # Normalise the spectrum
        if norm is None:
            data, norm_val = _normalise(data, domain, k_edge_last_x)
            # norm_val = 1
            # print("Not normalising")
        else:
            data, norm_val = _normalise(data, domain, k_edge_last_x, norm_val=norm)
            # norm_val = 1
            # print("Not normalising")

        return domain, data, norm_val


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
    plot_i_atoms,
):
    # Initialise user-defined arrays and variables
    # Broadening parameters
    # Start and end values of spectra
    start = begin
    stop = end

    # Energy of the lead peak
    first_peak = initial_peak

    # Set the start and end point of the linearly increasing broadening
    # change from the first and last ranges with respect to the leading
    # peak
    ewid_1 = first_peak + 5.0
    ewid_2 = first_peak + 15.0

    # System parameters ###
    # Chemical system
    # molecule = "graphene"
    # metal = "Cu"
    # element = "C"
    metal = surface
    element = excited_atom

    # Index range of the atom directories created by autoscript.py
    start_atom = index_start
    end_atom = index_end

    # Type of NEXAFS spectrum to output
    # 1 for total summed NEXAFS, 2 for angular, 3 for polarised, and
    # 4 for average polarised
    spectrum_type_opts = np.array(
        ["tot_summed", "angular", "polarised", "avg_polarised"]
    )
    n_type = np.where(spectrum_type_opts == spectrum_type)[0][0] + 1

    # The element index of the excited atom all the elements in the system
    # always the last element in the system, so if system contains H, C, Ag, C:exc
    # it will be 4
    atom = excited_nth_element

    # Create a list of all the atoms
    atom_ids = list(range(start_atom, end_atom + 1))

    # Set up a list of all the directories all the data is in C48/, C49/... C57/
    dirs = np.array([])
    for n in atom_ids:
        direc = f"{element}{str(n)}/"
        if os.path.isdir(direc):
            dirs = np.append(dirs, direc)

    # Deltas file for each molpdos run
    fname = f"/{molecule}_{metal}_{atom}_1_1_1_deltas.dat"

    # Get the length of the deltas file
    tmp_bands = np.loadtxt(f"{element}{str(atom_ids[0])}/t{theta[0]}_p{phi[0]}{fname}")
    tmp_bands_l = len(tmp_bands)

    # Create arrays with sizes of the system to use
    peaks = np.zeros([len(atom_ids), tmp_bands_l])
    bands = np.zeros([len(atom_ids), tmp_bands_l])

    # Plot spectrum for all theta and phi angles
    for t in theta:
        for p in phi:
            flat_peaks, flat_bands, peaks, bands = get_nexafs_data(
                t,
                p,
                dirs,
                fname,
                peaks,
                bands,
                n_type,
                molecule,
                metal,
                get_i_atom=plot_i_atoms,
            )

            print(f"Broadening delta peaks for theta={t} phi={p}...")
            print("Broadening total spectrum...")
            # x, y = broaden(
            x, y, norm_val = broaden(
                n_procs,
                start,
                stop,
                flat_peaks,
                flat_bands,
                omega_1,
                omega_2,
                mix_1,
                mix_2,
                ewid_1,
                ewid_2,
            )

            # Write out spectrum to a text file
            np.savetxt(
                f"{molecule}_{metal}_spectrum_t{t}_p{p}.txt", np.vstack((x, y)).T
            )

            if plot_i_atoms:
                print("Broadening individual atom spectra...")

                for i in range(len(dirs)):
                    print(f"Broadening atom {i+1}...")
                    # x, y = broaden(
                    x, y, _ = broaden(
                        n_procs,
                        start,
                        stop,
                        peaks[i, :],
                        bands[i, :],
                        omega_1,
                        omega_2,
                        mix_1,
                        mix_2,
                        ewid_1,
                        ewid_2,
                        with_tqdm=False,
                        norm=norm_val,
                    )

                    # Write out spectrum to a text file
                    np.savetxt(
                        f"{dirs[i]}t{t}_p{p}/{molecule}_{metal}_spectrum_t{t}_p{p}.txt",
                        np.vstack((x, y)).T,
                    )

                print("Finished broadening individual atom spectra...")


@click.command()
@click.option(
    "-n",
    "--nprocs",
    default=1,
    type=bool,
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
    "-s",
    "--spectrum_type",
    default="avg_polarised",
    type=click.Choice(["tot_summed", "angular", "polarised", "avg_polarised"]),
    help="type of spectrum to plot",
)
@click.option("-p", "--phi", default=np.array([60]), type=np.array, help="phi angles")
@click.option(
    "-t",
    "--theta",
    default=np.array([0, 25, 53, 90]),
    type=np.array,
    help="theta angles",
)
@click.option(
    "-x",
    "--excited_nth_element",
    default=3,
    type=int,
    help="the element index of the excited atom of all the elements in the system - this will always be the last element in the system, so if system contains C, Cu, C:exc it will be 3",
)
@click.option(
    "-p",
    "--plot_i_atoms",
    is_flag=True,
    default=True,
    type=bool,
    help="plot individual atom spectra",
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
    plot_i_atoms,
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
        plot_i_atoms,
    )


if __name__ == "__main__":
    nexafs()
