#!/usr/bin/env python
# coding: utf-8
"""Process global UM output by interpolating selected fields to common grid."""

# Commonly used standard library tools
import argparse
from pathlib import Path
from time import time
import warnings

# My packages and local scripts
from aeolus.core import Run

from commons import (
    GLM_MODEL_TIMESTEP,
    GLM_START_DAY,
)
import mypaths
from proc_um_output import get_filename_list, process_cubes
from utils import create_logger


# Global definitions and styles
warnings.filterwarnings("ignore")
SCRIPT = Path(__file__).name


def parse_args(args=None):
    """Argument parser."""
    ap = argparse.ArgumentParser(
        SCRIPT,
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    ap.add_argument(
        "-p", "--planet", type=str, required=True, help="Planet configuration"
    )
    ap.add_argument("-r", "--run", type=str, required=True, help="Run key")
    ap.add_argument(
        "-s",
        "--startday",
        type=int,
        default=GLM_START_DAY,
        help="Load files with timestamp >= this",
    )

    return ap.parse_args(args)


def main(args=None):
    """Main entry point of the script."""
    # Parse command-line arguments
    args = parse_args(args)
    planet = args.planet
    run_key = args.run

    label = f"{planet}_{run_key}"
    L.info(f"label = {label}")

    # Create a subdirectory for processed data
    outdir = mypaths.sadir / label / "_processed"
    outdir.mkdir(parents=True, exist_ok=True)

    # Make a list of files matching the file mask and the start day threshold
    fnames = get_filename_list(mypaths.sadir / label, ts_start=args.startday)
    L.debug(f"fnames = {fnames[0]} ... {fnames[-1]}")

    # Initialise a `Run` by loading data from the selected files
    run = Run(files=fnames, name=label, planet=planet, timestep=GLM_MODEL_TIMESTEP)

    # Regrid & interpolate data
    run.proc_data(process_cubes, timestep=run.timestep)

    # Write the data to a netCDF file
    fname_out = outdir / f"{run.name}.nc"
    run.to_netcdf(fname_out)
    L.success(f"Saved to {fname_out}")


if __name__ == "__main__":
    t0 = time()
    L = create_logger(Path(__file__))
    main()
    L.info(f"Execution time: {time() - t0:.1f}s")
