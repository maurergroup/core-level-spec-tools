#!/usr/bin/env python3

import os

import numpy as np
from tqdm import tqdm


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


def broaden(
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
        for i in tqdm(range(len(domain))):
            data[i] = np.sum(
                _schmid_pseudo_voigt(domain[i], mixing, dirac_peaks, sigma) * coeffs
            )

    else:
        for i in range(len(domain)):
            data[i] = np.sum(
                _schmid_pseudo_voigt(domain[i], mixing, dirac_peaks, sigma) * coeffs
            )

    # Normalise the spectrum
    if norm is None:
        data, norm_val = _normalise(data, domain, k_edge_last_x)
    else:
        data, norm_val = _normalise(data, domain, k_edge_last_x, norm_val=norm)

    return domain, data, norm_val


def main():
    # Initialise user-defined arrays and variables
    # Broadening parameters
    # Start and end values of spectra
    start = 287.0
    stop = 304.0

    # Broadening parameters for the first and last ranges
    omega_1 = 0.75
    omega_2 = 2.0

    # Energy of the lead peak
    first_peak = 290.0

    # Set the start and end point of the linearly increasing broadening
    # change from the first and last ranges with respect to the leading
    # peak
    ewid_1 = first_peak + 5.0
    ewid_2 = first_peak + 15.0

    # Set the Gaussian/Lorentzian mixing ratio for the ranges
    mix_1 = 0.2
    mix_2 = 0.8

    # System parameters ###
    # Chemical system
    molecule = "graphene"
    metal = "Cu"
    element = "C"

    # Index range of the atom directories created by autoscript.py
    start_atom = 234
    end_atom = 454

    # Type of NEXAFS spectrum to output
    # 1 for total summed NEXAFS, 2 for angular, 3 for polarised, and
    # 4 for average polarised
    n_type = 4

    # The theta and phi angles simulated
    theta = np.array(["00", "25", "53", "90"])
    phi = np.array(["60"])

    # The element index of the excited atom all the elements in the system
    # always the last element in the system, so if system contains H, C, Ag, C:exc
    # it will be 4
    atom = "3"

    # Plot individual atom spectra
    plot_i_atoms = True

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


if __name__ == "__main__":
    main()
