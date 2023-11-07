#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np


def plot(angle):
    x, y = np.loadtxt(f"graphene_Cu_spectrum_{angle}.txt", usecols=(0, 1), unpack=True)
    # layers = np.loadtxt("layers.txt")

    # dirac_peak_colours = ["#FF2F92", "#00F900", "#0096FF", "#9437FF", "#000000"]
    # layer_names = [
    #     "defect atom",
    #     "1 nearest neighbour",
    #     "2 nearest neighbours",
    #     "3 nearest neighbours",
    #     "4 nearest neighbours",
    # ]
    # first_call = [True, True, True, True, True]

    # # Plot the individual binding energies
    # for peak in x:
    #     # Find the layer the peak is in
    #     layer = np.where(layers == float(peak.split()[1]))[1][0]

    #     # Include the peak in the lengend if first call
    #     if first_call[layer] is True:
    #         plt.axvline(
    #             x=float(peak.split()[0]),
    #             c=dirac_peak_colours[layer],
    #             ymax=0.25,
    #             label=layer_names[layer],
    #         )
    #         first_call[layer] = False
    #     else:
    #         plt.axvline(
    #             x=float(peak.split()[0]), c=dirac_peak_colours[layer], ymax=0.25
    #         )

    for peak in x:
        plt.axvline(x=float(peak), ymax=0.25)

    plt.plot(x, y)
    plt.savefig(f"spectrum_{angle}.png")


def main():
    deltas = ["t00_p60", "t25_p60", "t53_p60", "t90_p60"]

    for i in deltas:
        plot(i)


if __name__ == "__main__":
    main()
