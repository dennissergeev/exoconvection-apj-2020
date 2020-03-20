#!/usr/bin/env python
# coding: utf-8
"""Same as `proc_global_um_output.py`, but for multiple folders within the ensemble."""

# Commonly used standard library tools
import configparser
from pathlib import Path
from time import time
import warnings

# My packages and local scripts
from aeolus.core import Run

from commons import GLM_MODEL_TIMESTEP
import mypaths
from proc_global_um_output import get_filename_list, parse_args, process_cubes
from utils import create_logger, pprint_dict


# Global definitions and styles
warnings.filterwarnings("ignore")
SCRIPT = Path(__file__).name


def main(args=None):
    """Main entry point of the script."""
    # Parse command-line arguments
    args = parse_args(args)
    planet = args.planet
    topdir = mypaths.sadir / f"{planet}_grcs_ensemble"
    labels = ["base"]
    labels += [
        i.with_suffix("").name.replace("rose-app-", "")
        for i in sorted(Path(str(topdir) + "_conf").glob("rose-app-*"))
    ]
    L.debug(f"labels = {labels}")
    for label in labels:
        L.info(f"label = {label}")
        if label == "base":
            config_str = ""
        else:
            config = configparser.ConfigParser()
            config.read(Path(str(topdir) + "_conf") / f"rose-app-{label}.conf")
            config_str = pprint_dict(dict(config["namelist:run_convection"]))
        # Make a list of files matching the file mask and the start day threshold
        fnames = get_filename_list(topdir / label, ts_start=args.startday)
        L.debug(f"fnames = {fnames[0]} ... {fnames[-1]}")

        # Create a subdirectory for processed data
        outdir = topdir / label / "_processed"
        outdir.mkdir(parents=True, exist_ok=True)

        # Initialise a `Run` by loading data from the selected files
        run = Run(
            files=fnames,
            description=config_str,
            name=label,
            planet=planet,
            timestep=GLM_MODEL_TIMESTEP,
        )

        # Regrid & interpolate data
        run.proc_data(process_cubes, timestep=run.timestep)

        # Write the data to a netCDF file
        fname_out = outdir / f"{planet}_{run.name}.nc"
        run.to_netcdf(fname_out)
        L.success(f"Saved to {fname_out}")


if __name__ == "__main__":
    t0 = time()
    L = create_logger(Path(__file__))
    main()
    L.info(f"Execution time: {time() - t0:.1f}s")
