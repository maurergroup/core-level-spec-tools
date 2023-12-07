#!/usr/bin/env python3

import os

import click


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
def molpdos():
    """
    Wrapper for execute_molpdos.sh
    """
    # Check that MolPDOS exists and is executable

    # Check that the individual atomic files exist

    return os.system(
        f"{os.path.dirname(os.path.realpath(__file__))}/execute_molpdos.sh"
    )
