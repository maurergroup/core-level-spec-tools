#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["figure.figsize"] = (8, 5)


def _plot(angle, path_1, label, path_2=None, colour=None):
    if path_2 is not None:
        x, y = np.loadtxt(
            f"{path_2}/graphene_Cu_spectrum_{angle}.txt", usecols=(0, 1), unpack=True
        )

        plt.plot(x, y, c=colour[1], label=label[1])

        colour = colour[0]
        label = label[0]

    x, y = np.loadtxt(
        f"{path_1}/graphene_Cu_spectrum_{angle}.txt", usecols=(0, 1), unpack=True
    )

    plt.plot(x, y, c=colour, label=label)

    return x, y


def multi_element(deltas, cmap):
    start_atom = 234
    end_atom = 403
    element = "C"

    # Create a list of all the atoms
    atom_ids = list(range(start_atom, end_atom + 1))

    # Set up a list of all the directories all the data is in C48/, C49/... C57/
    dirs = np.array([])
    for n in atom_ids:
        direc = f"{element}{str(n)}/"
        if os.path.isdir(direc):
            dirs = np.append(dirs, direc)

    for i in deltas:
        # for j, c in zip(dirs, np.linspace(0.2, 0.6, len(dirs))):
        for j, c in zip(dirs, np.linspace(0.1, 1, len(dirs))):
            _plot(i, f"{j}{i}", j[:-1], cmap(c))
            # plot(i, f"{j}{i}", j[:-1])

        # x, y = plot(i, "./", "Total Spectrum", "blueviolet")
        x, _ = _plot(i, "./", "Total Spectrum", "black")
        plt.xticks(np.arange(min(x), max(x) + 1, 2))
        plt.ylabel("Normalised Intensity")
        plt.xlabel("Energy (eV)")
        plt.tight_layout()
        plt.legend()

        plt.savefig(f"spectrum_{i}.png")
        print(f"Plotted {i}")

        plt.figure().clear()


def multi_angle(deltas, cmap):
    for delta, c in zip(deltas, np.linspace(0.1, 1, len(deltas))):
        x, _ = _plot(delta, "./", rf"$\theta = {delta[1:3]}$", cmap(c))

    plt.xticks(np.arange(min(x), max(x) + 1, 2))
    plt.ylabel("Intensity")
    plt.xlabel("Energy (eV)")
    plt.tight_layout()
    plt.legend()

    plt.savefig(f"spectrum_angles.png")
    print(f"Plotted angles")


def compare(deltas):
    prist_path = "./pristine/setup/NEXAFS"
    sw_path = "./sw/setup/2_NEXAFS"

    for delta in deltas:
        x, _ = _plot(
            delta,
            prist_path,
            ["Pristine", "Stone-Wales"],
            path_2=sw_path,
            colour=["darkviolet", "darkorange"],
        )

        plt.xticks(np.arange(min(x), max(x) + 1, 2))
        plt.ylabel("Intensity")
        plt.xlabel("Energy (eV)")
        plt.tight_layout()
        plt.legend()

        plt.savefig(f"spectrum_{delta}.png")
        print(f"Plotted {delta}")

        plt.figure().clear()


def main():
    # NEXAFS angles
    deltas = ["t00_p60", "t25_p60", "t53_p60", "t90_p60"]

    # Use a linear segmented colormap
    cmap = plt.get_cmap("plasma").reversed()

    # multi_element(deltas, cmap)
    # multi_angle(deltas, cmap)
    compare(deltas)


if __name__ == "__main__":
    main()
