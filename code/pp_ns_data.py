#!/usr/bin/env python
# coding: utf-8
"""Process NS UM output by interpolating selected fields to common grid."""

# Commonly used standard library tools
import argparse
from datetime import timedelta
from pathlib import Path
from time import time
import warnings

# My packages and local scripts
from aeolus.core import Run

from commons import (
    DT_FMT,
    NS_CYCLE_TIMES,
    NS_MODEL_TYPES,
)
import mypaths
from proc_um_output import process_cubes
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

    return ap.parse_args(args)


def main(args=None):
    """Main entry point of the script."""
    # Parse command-line arguments
    args = parse_args(args)
    planet = args.planet
    L.info(f"planet = {planet}")
    run_key = args.run
    L.info(f"run_key = {run_key}")
    subdir = f"{planet}_{run_key}"

    label = f"{planet}_{run_key}"
    L.info(f"label = {label}")
    for model_type, model_specs in NS_MODEL_TYPES.items():
        L.info(f"model_type = {model_type}")
        label = f"{planet}_{run_key}_{model_type}"
        for i in range(NS_CYCLE_TIMES[subdir]["ndays"]):
            _cycle = (NS_CYCLE_TIMES[subdir]["start"] + timedelta(days=i)).strftime(
                DT_FMT
            )
            L.info(f"_cycle = {_cycle}")
            fpath = mypaths.nsdir / subdir / _cycle / model_specs["path"]
            L.info(f"fpath = {fpath}")
            # Create a subdirectory for processed data
            outdir = fpath.parent / "_processed"
            outdir.mkdir(parents=True, exist_ok=True)

            # Initialise a `Run` by loading data from the selected files
            run = Run(
                files=fpath,
                name=label,
                planet=planet,
                model_type=model_type,
                timestep=model_specs["timestep"],
            )

            # Regrid & interpolate data
            run.proc_data(process_cubes, timestep=run.timestep)

            # Write the data to a netCDF file
            fname_out = outdir / f"{run.name}_{_cycle}.nc"
            run.to_netcdf(fname_out)
            L.success(f"Saved to {fname_out}")


if __name__ == "__main__":
    t0 = time()
    L = create_logger(Path(__file__))
    main()
    L.info(f"Execution time: {time() - t0:.1f}s")
