#!/usr/bin/env python
# coding: utf-8
"""Process global UM output by interpolating selected fields to common grid."""

# Commonly used standard library tools
import argparse
from pathlib import Path
import re
from time import time
import warnings

# Scientific stack
import iris

# My packages and local scripts
from aeolus.coord_utils import (
    regrid_3d,
    replace_z_coord,
    ensure_bounds,
    UM_HGT,
)
from aeolus.core import Run
from aeolus.grid import roll_cube_pm180
from aeolus.subset import CM_MEAN_CONSTR

import mypaths
from utils import create_logger


# Global definitions and styles
warnings.filterwarnings("ignore")
SCRIPT = Path(__file__).name

# Common paths
run_keys = ["base"]  # use only base run

MODEL_TIMESTEP = 86400 / 72  # in seconds
RUNID = r"umglaa"  # file prefix
FILE_REGEX = RUNID + r".p[b,c,d,e]{1}[0]{6}(?P<timestamp>[0-9]{2,4})_00"
DEFAULT_START_DAY = 1440  # by default, use only the last 360 days of a 5 year run

# Selected variables
SINGLE_LEVEL_VARS = [
    "surface_temperature",
    "toa_incoming_shortwave_flux",
    "toa_outgoing_longwave_flux",
    "toa_outgoing_longwave_flux_assuming_clear_sky",
    "toa_outgoing_shortwave_flux",
    "toa_outgoing_shortwave_flux_assuming_clear_sky",
    "convective_rainfall_flux",
    "convective_snowfall_flux",
    "high_type_cloud_area_fraction",
    "low_type_cloud_area_fraction",
    "medium_type_cloud_area_fraction",
    "precipitation_flux",
    "stratiform_rainfall_flux",
    "stratiform_snowfall_flux",
]
MULTI_LEVEL_VARS = [
    "air_potential_temperature",
    "air_pressure",
    "specific_humidity",
    "mass_fraction_of_cloud_liquid_water_in_air",
    "mass_fraction_of_cloud_ice_in_air",
    "upward_air_velocity",
]
HORIZ_WINDS_CONSTR = ["x_wind", "y_wind"]
INCR_CONSTR = iris.AttributeConstraint(STASH=lambda x: x.item in [181, 182, 233])


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
        default=DEFAULT_START_DAY,
        help="Load files with timestamp >= this",
    )

    return ap.parse_args(args)


def process_cubes(cubelist, timestep=MODEL_TIMESTEP, ref_cube_constr="specific_humidity"):
    """Post-process data for easier analysis."""
    cubes = iris.cube.CubeList()

    # First, extract all multi-level fields (30-day averages)
    cubes += cubelist.extract(MULTI_LEVEL_VARS).extract(CM_MEAN_CONSTR)

    # Horizontal wind components
    for cube in cubelist.extract(HORIZ_WINDS_CONSTR):
        cube.units = "m s^-1"
        cubes.append(cube)

    # Increments
    for cube in cubelist.extract(INCR_CONSTR):
        if cube.attributes["STASH"].item in [181, 233]:
            incr_unit = "K"
        else:
            incr_unit = "kg kg^-1"
        if cube.attributes["STASH"].item == 233:
            cube.units = f"{incr_unit} s^-1"
        else:
            cube.units = f"{1/timestep} {incr_unit} s^-1"
            cube.convert_units(f"{incr_unit} s^-1")
        cubes.append(cube)

    # Interpolation & regridding to common grid
    ref_cube = cubes.extract_strict(ref_cube_constr)
    ref_cube = replace_z_coord(ref_cube, UM_HGT)

    # Interpolate to common levels
    cubes = iris.cube.CubeList(
        [regrid_3d(replace_z_coord(cube, UM_HGT), ref_cube, UM_HGT) for cube in cubes]
    )

    # Add all single-level cubes
    cubes += cubelist.extract(SINGLE_LEVEL_VARS).extract(CM_MEAN_CONSTR)

    # Roll cubes to +/- 180 degrees in longitude for easier analysis
    r_cubes = iris.cube.CubeList()
    for cube in cubes:
        r_c = roll_cube_pm180(cube)
        ensure_bounds(r_c)
        r_cubes.append(r_c)

    return r_cubes


def main(args=None):
    """Main entry point of the script."""
    # Parse command-line arguments
    args = parse_args(args)
    planet = args.planet
    run_key = args.run

    label = f"{planet}_{run_key}"

    # Create a subdirectory for processed data
    outdir = mypaths.sadir / label / "_processed"
    outdir.mkdir(parents=True, exist_ok=True)

    # Make a list of files matching the file mask and the start day threshold
    fnames = []
    for fpath in sorted((mypaths.sadir / label).glob(RUNID + "*")):
        match = re.match(FILE_REGEX, fpath.name,)
        if match:
            if int(match["timestamp"]) >= args.startday:  #
                fnames.append(fpath)
    # Initialise a `Run` by loading data from the selected files
    run = Run(files=fnames, name=label, planet=planet, timestep=MODEL_TIMESTEP)

    # Regrid & interpolate data
    run.proc_data(process_cubes, timestep=run.timestep)

    # Remove planet_conf attribute before saving
    cubes_out = iris.cube.CubeList()
    for cube in run.proc:
        try:
            cube.attributes.pop("planet_conf")
        except KeyError:
            pass
        cubes_out.append(cube)

    # Write the data to a netCDF file
    fname_out = outdir / f"{run.name}.nc"
    iris.save(cubes_out, str(fname_out))
    L.success(f"Saved to {fname_out}")


if __name__ == "__main__":
    t0 = time()
    L = create_logger(Path(__file__))
    main()
    L.info(f"Execution time: {time() - t0:.1f}s")
