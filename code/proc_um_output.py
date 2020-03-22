#!/usr/bin/env python
# coding: utf-8
"""Functions for post-processing UM output."""

# Commonly used standard library tools
import re
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
from aeolus.grid import roll_cube_pm180
from aeolus.subset import CM_MEAN_CONSTR

from commons import (
    GLM_FILE_REGEX,
    GLM_MODEL_TIMESTEP,
    GLM_RUNID,
    INCR_CONSTR,
    HORIZ_WINDS_NAMES,
    HORIZ_WINDS_STASH,
    MULTI_LEVEL_VARS,
    SINGLE_LEVEL_VARS,
)


def process_cubes(
    cubelist,
    timestep=GLM_MODEL_TIMESTEP,
    ref_cube_constr="specific_humidity",
    extract_mean=True,
):
    """Post-process data for easier analysis."""
    cm_constr = iris.Constraint()
    if extract_mean:
        cm_constr &= CM_MEAN_CONSTR

    cubes = iris.cube.CubeList()

    # First, extract all multi-level fields (30-day averages)
    cubes += cubelist.extract(MULTI_LEVEL_VARS).extract(cm_constr)

    # Horizontal wind components
    winds = cubelist.extract(HORIZ_WINDS_NAMES)
    if len(winds) < 2:
        winds = cubelist.extract(
            iris.AttributeConstraint(STASH=lambda x: x in HORIZ_WINDS_STASH)
        )
    if len(winds) == 2:
        for (cube, name) in zip(winds, HORIZ_WINDS_NAMES):
            cube.units = "m s^-1"
            cube.rename(name)
            cubes.append(cube)
    else:
        warnings.warn("Unable to extract unique cubes for horizontal winds")

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
    cubes += cubelist.extract(SINGLE_LEVEL_VARS).extract(cm_constr)

    # Roll cubes to +/- 180 degrees in longitude for easier analysis
    r_cubes = iris.cube.CubeList()
    for cube in cubes:
        r_c = roll_cube_pm180(cube)
        ensure_bounds(r_c)
        r_cubes.append(r_c)

    return r_cubes


def get_filename_list(
    path_to_dir,
    glob_pattern=f"{GLM_RUNID}*",
    ts_start=0,
    regex=GLM_FILE_REGEX,
    regex_key="timestamp",
):
    """Get a list of files with timestamps greater or equal than start in a directory."""
    fnames = []
    for fpath in sorted(path_to_dir.glob(glob_pattern)):
        match = re.match(regex, fpath.name)
        if match:
            if int(match[regex_key]) >= ts_start:
                fnames.append(fpath)
    return fnames
