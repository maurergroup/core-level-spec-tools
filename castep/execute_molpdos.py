#!/usr/bin/env python3

import os


def run():
    """
    Wrapper for execute_molpdos.sh
    """
    return os.system(
        f"{os.path.dirname(os.path.realpath(__file__))}/execute_molpdos.sh"
    )
