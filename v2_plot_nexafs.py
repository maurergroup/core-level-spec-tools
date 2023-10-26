#!/usr/bin/env python3

import os

import numpy as np
from tqdm import tqdm


def get_nexafs_data(
    theta, phi, dirs, fname, peaks, bands, n_type, molecule, metal
) -> tuple[np.ndarray, np.ndarray]:
    """Parse the NEXAFS data for all theta and phi"""

    for i in range(len(dirs)):
        nexafs = np.loadtxt(f"{dirs[i]}t{theta}_p{phi}{fname}")
        peaks[i, :] = nexafs[:, 0]
        bands[i, :] = nexafs[:, n_type]

    print(f"Finished parsing NEXAFS data for theta={theta} phi={phi}")

    flat_peaks = peaks.flatten()
    flat_bands = bands.flatten()

    # Write out all of the data into a delta peaks file
    np.savetxt(
        f"{molecule}_{metal}_deltas_t{theta}_p{phi}.txt",
        np.vstack((flat_peaks, flat_bands)).T,
        header="# <x in eV> Intensity",
    )

    print(f"Finished writing out delta peaks file for theta={theta} phi={phi}")

    return flat_peaks, flat_bands


def _schmid_pseudo_voigt(
    domain, m, E, omega, asymmetry=False, a=None, b=None
) -> np.ndarray:
    """
    Apply broadening scheme for XPS spectra
    https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/sia.5521

    domain = linspace of x range and bin width
    A = intensity
    m = Gaussian-Lorentzian mixing parameter
    E = line centre (aka dirac peak)
    omega = full width at half maximum (omega > 0)
    asymmetry = True or False
    a = asymmetry parameter
    b = asymmetry translation parameter
    """

    # if asymmetry is False:
    # A has been omitted for performance as it will almost certainly only be set to 1
    return (1 - m) * np.sqrt((4 * np.log(2)) / (np.pi * omega**2)) * np.exp(
        -(4 * np.log(2) / omega**2) * (domain - E) ** 2
    ) + m * (1 / (2 * np.pi)) * (omega / ((omega / 2) ** 2 + (domain - E) ** 2))

    # else:
    #     omega_as = 2 * omega / (1 + np.exp(-a * domain - b))

    #     return A * (1 - m) * np.sqrt(
    #         (4 * np.log(2))
    #         / (np.pi * ((2 * omega_as) / (1 + np.exp(-a * ((domain - E) - b)))) ** 2)
    #     ) * np.exp(
    #         -(
    #             4
    #             * np.log(2)
    #             / (2 * omega_as / (1 + np.exp(-a * ((domain - E) - b)))) ** 2
    #         )
    #         * (domain - E) ** 2
    #     ) + A * m * (
    #         1 / (2 * np.pi)
    #     ) * (
    #         (2 * omega_as / (1 + np.exp(-a * ((domain - E) - b))))
    #         / (
    #             (((2 * omega_as) / (1 + np.exp(-a * ((domain - E) - b)))) / 2) ** 2
    #             + (domain - E) ** 2
    #         )
    #     )


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
) -> tuple[np.ndarray, np.ndarray]:
    """
    Broaden dirac delta peaks

    start = beginning of x range
    stop = end of x range
    dirac_peaks = list of dirac delta peaks
    eta = full width at half maximum (eta > 0)
    """
    domain = np.arange(start, stop, bin_width)
    data = np.zeros([len(domain)])

    if coeffs is None:
        coeffs = np.zeros(len(dirac_peaks))

    # Find the peaks in the different broadening regions
    mic_1_trans = np.where(dirac_peaks <= ewid_1)
    mic_2_trans = np.where((dirac_peaks > ewid_1) & (dirac_peaks <= ewid_2))
    m1t_len = len(mic_1_trans)
    m2t_len = len(mic_2_trans)

    sigma = np.zeros((len(dirac_peaks)))
    mixing = np.zeros((len(dirac_peaks)))

    for i in range(m1t_len):
        sigma[i] = omega_1
        mixing[i] = mix_1

    for i in np.arange(m1t_len, m1t_len + m2t_len):
        sigma[i] = omega_2
        mixing[i] = mix_2

    sigma_precalc = omega_1 + ((omega_2 - omega_1) / (ewid_2 - ewid_1))
    mix_precalc = mix_1 + (mix_2 - mix_1) / (ewid_2 - ewid_1)

    for i in np.arange(m1t_len + m2t_len, len(dirac_peaks)):
        sigma[i] = sigma_precalc * (dirac_peaks[i] - ewid_1)
        mixing[i] = mix_precalc * (dirac_peaks[i] - ewid_1)

    for i in tqdm(range(len(domain))):
        V = _schmid_pseudo_voigt(domain[i], mixing, dirac_peaks, sigma) * coeffs
        data[i] = np.sum(V)

    return domain, data


def main():
    # Initialise user-defined arrays and variables
    # Broadening parameters
    # Start and end values of spectra
    start = 285.0
    stop = 300.0

    # Broadening parameters for the first and last ranges
    omega_1 = 0.75
    omega_2 = 2.0

    # Energy of the lead peak
    first_peak = 289.0

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
    start_atom = 248
    end_atom = 477

    # Type of NEXAFS spectrum to output
    # 1 for total summed NEXAFS, 2 for angular, 3 for polarised, and
    # 4 for average polarised
    n_type = 1

    # The theta and phi angles simulated
    theta = np.array(["00", "25", "53", "90"])
    phi = np.array(["60"])

    # The element index of the excited atom all the elements in the system
    # always the last element in the system, so if system contains H, C, Ag, C:exc
    # it will be 4
    atom = "3"

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

    # Create arrays with sizes of the system to use
    peaks = np.zeros([len(atom_ids), len(tmp_bands)])
    bands = peaks

    # Plot spectrum for all theta and phi angles
    for t in theta:
        for p in phi:
            peaks, bands = get_nexafs_data(
                t, p, dirs, fname, peaks, bands, n_type, molecule, metal
            )

            print(f"Broadening delta peaks for theta={t} phi={p}...")

            x, y = broaden(
                start,
                stop,
                peaks,
                bands,
                omega_1,
                omega_2,
                mix_1,
                mix_2,
                ewid_1,
                ewid_2,
            )

            print(f"Finished broadening delta peaks for theta={t} phi={p}")

            # Write out spectrum to a text file
            np.savetxt(f"{molecule}_{metal}_spectrum_t{t}_p{p}.txt", np.array(x, y))


if __name__ == "__main__":
    main()
