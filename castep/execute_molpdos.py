#!/usr/bin/env python3

import os
import re
from shutil import which

import click


def _check_executable(executable):
    """
    Check if an executable is available.

    Parameters
    ----------
        executable : str
            Name of the executable to check
    """

    if which(executable) is None:
        raise FileNotFoundError(
            f"{executable} not found. Please make sure it is installed and available "
            "in $PATH."
        )


def _check_atom_dirs(element):
    """
    Check that the individual atomic files exist for the given element.

    Parameters
    ----------
        element : str
            Element symbol of the excited atom
    """

    # Get a list of directories in the current working directory
    dirs = [d for d in os.listdir(".") if os.path.isdir(d)]

    # element followed by 1-4 digits
    # I will be amazed if anyone can run NEXAFS with >9999 atoms
    regexs = []
    for dir in dirs:
        regexs.append(re.search(element + "[0-9]{1,4}$", dir))

    # Error if no directories match
    if len(regexs) == 0:
        raise FileNotFoundError(f"No directories found for {element}.")


def main(phi, theta, molecule, surface, excited_atom):
    """
    Wrapper for execute_molpdos.sh

    Parameters
    ----------
        phi : list[str]
            phi angles
        theta : list[str]
            theta angles
        molecule : str
            adsorbate name
        surface : str
            element that comprises the surface
        excited_atom : str
            element symbol of the excited atom
    """

    # Check that the required executables are available
    _check_executable("molpdos")

    # Check that the required directories exist
    _check_atom_dirs(excited_atom)

    # Convert phi and theta args into bash-compatible versions
    phi_b = " ".join(phi)
    theta_b = " ".join(theta)

    return os.system(
        f"{os.path.dirname(os.path.realpath(__file__))}/execute_molpdos.sh "
        f"{phi_b} {theta_b} {molecule} {surface} {excited_atom}"
    )


@click.command()
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
@click.option("-m", "--molecule", required=True, type=str, help="adsorbate name")
@click.option("-s", "--surface", required=True, type=str, help="surface element")
@click.option(
    "-a",
    "--excited_atom",
    required=True,
    type=str,
    help="element symbol of the excited atom",
)
def molpdos(phi, theta, molecule, surface, excited_atom):
    """
    Calculate MolPDOS for excited atoms in a CASTEP NEXAFS calculation.

    Copyright \u00A9 2023-2024, Dylan Morgan dylan.morgan@warwick.ac.uk
    """

    main(phi, theta, molecule, surface, excited_atom)
