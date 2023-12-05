#!/usr/bin/env python3

import os


def run():
    """
    Wrapper for execute_molpdos.sh
    """
    return os.system(f"{os.getcwd()}/castep/execute_molpdos.sh")
