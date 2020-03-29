#!/usr/bin/env python
# coding: utf-8

"""Extracting mean vertical profiles from the nesting suite."""

# Commonly used standard library tools
from datetime import timedelta
from pathlib import Path
from time import time
import warnings

# Scientific stack
import iris

import xarray as xr

# My packages and local scripts
from aeolus.calc import spatial
from aeolus.core import Run
from aeolus.subset import DIM_CONSTR_ZYX

from commons import (
    DT_FMT,
    NS_CYCLE_TIMES,
    NS_MODEL_TYPES,
    NS_RUN_ALIASES,
    PLANET_ALIASES,
    SS_REGION,
)
from gl_diag import (
    calc_derived_cubes,
)
import mypaths
from utils import create_logger


warnings.filterwarnings("ignore")
L = create_logger(Path(__file__))


@L.catch
def main():
    """Main function."""
    # Create a subdirectory for processed data
    outdir = mypaths.nsdir / "_processed"
    outdir.mkdir(parents=True, exist_ok=True)

    for planet in PLANET_ALIASES.keys():
        L.info(f"planet = {planet}")
        for run_key in NS_RUN_ALIASES.keys():
            L.info(f"run_key = {run_key}")
            subdir = f"{planet}_{run_key}"
            for model_type, model_specs in NS_MODEL_TYPES.items():
                L.info(f"model_type = {model_type}")
                label = f"{planet}_{run_key}_{model_type}"
                tmp_cl = iris.cube.CubeList()
                for i in range(NS_CYCLE_TIMES[subdir]["ndays"]):
                    _cycle = (
                        NS_CYCLE_TIMES[subdir]["start"] + timedelta(days=i)
                    ).strftime(DT_FMT)
                    L.info(f"_cycle = {_cycle}")
                    fpath = (
                        mypaths.nsdir
                        / subdir
                        / _cycle
                        / model_specs["path"].parent
                        / "_processed"
                        / f"{label}_{_cycle}.nc"
                    )
                    L.info(f"fpath = {fpath}")
                    # Initialise a `Run` by loading processed data
                    run = Run(
                        files=fpath,
                        name=label,
                        planet=planet,
                        model_type=model_type,
                        timestep=model_specs["timestep"],
                        processed=True,
                    )
                    # Derive additional fields
                    run.add_data(calc_derived_cubes)

                    for cube in run.proc.extract(DIM_CONSTR_ZYX):
                        L.info(f"cube = {cube.name()}")
                        ave_cube = spatial(cube.extract(SS_REGION.constraint), "mean")
                        tmp_cl.append(ave_cube)

                for cube in tmp_cl:
                    try:
                        cube.attributes.pop("planet_conf")
                    except KeyError:
                        pass
                tmp_cl = tmp_cl.merge()
                L.debug(f"merged cubelist = {tmp_cl}")
                _dict = {}
                for cube in tmp_cl:
                    # cube.attributes = {
                    #     k: v for k, v in cube.attributes.items() if k != "planet_conf"
                    # }
                    _dict[cube.name()] = xr.DataArray.from_iris(cube)
                ds = xr.merge([_dict], compat="override")
                L.debug(f"dataset = {ds}")
                ds.attrs.update(
                    {
                        "region": SS_REGION.to_str(),
                        "model_type": run.model_type,
                        "run_key": run_key,
                        "planet": planet,
                    },
                )
                fname_out = outdir / f"ns_area_mean_vprof_{run.name}.nc"
                ds.to_netcdf(fname_out)
                L.success(f"Saved to {fname_out}")


if __name__ == "__main__":
    t0 = time()
    main()
    L.info(f"Execution time: {time() - t0:.1f}s")
