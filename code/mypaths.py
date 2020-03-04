# -*- coding: utf-8 -*-
"""Paths to data."""
from pathlib import Path

# Top-level directory containing code and data (one level up)
topdir = Path(__file__).absolute().parent.parent

# Modelling results
datadir = topdir / "modelling" / "um" / "results"

nsdir = datadir / "ns"  # nesting suites directory
sadir = datadir / "sa"  # standalone suites directory

# Plotting output
plotdir = topdir / "plots"


def lsdir(path):
    return sorted(path.glob("*"))
