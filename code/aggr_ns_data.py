#!/usr/bin/env python
# coding: utf-8
"""Aggregated metrics of NS data."""

# Standard library tools
from datetime import timedelta
from pathlib import Path
from time import time
import warnings

# Scientific stack
import iris

import xarray as xr

# My packages and local scripts
from aeolus.coord_utils import UM_TIME
from aeolus.core import Run
from aeolus.exceptions import MissingCubeError

from commons import (
    DT_FMT,
    NS_CYCLE_TIMES,
    NS_MODEL_TYPES,
    NS_RUN_ALIASES,
    PLANET_ALIASES,
    SS_REGION,
)
from gl_diag import (
    DIAGS,
    ONLY_GLOBAL,
    ONLY_LAM,
    calc_derived_cubes,
)
import mypaths
from utils import create_logger


warnings.filterwarnings("ignore")
L = create_logger(Path(__file__))


@L.catch
def main(args=None):
    """Main function."""

    def _do_calc(_calc, cl, cl_out):
        try:
            scalar_coords = cl[0].coords(UM_TIME)
            attrs = {
                k: v
                for k, v in cl[0].attributes.items()
                if k not in ["STASH", "planet_conf"]
            }
            cube = _calc(cl)
            cube.attributes = attrs
            for sc_c in scalar_coords:
                try:
                    cube.add_aux_coord(sc_c)
                except ValueError:
                    pass
            cl_out.append(cube)
        except (ValueError, MissingCubeError, iris.exceptions.ConstraintMismatchError) as e:
            L.debug(f"Caught exception:\n{e}\n")
            pass

    # Create a subdirectory for processed data
    outdir = mypaths.nsdir / "_processed"
    outdir.mkdir(parents=True, exist_ok=True)

    lam_diags = {k: v for k, v in DIAGS.items() if k not in ONLY_GLOBAL}

    for planet in PLANET_ALIASES.keys():
        L.info(f"planet = {planet}")
        for run_key in NS_RUN_ALIASES.keys():
            L.info(f"run_key = {run_key}")
            subdir = f"{planet}_{run_key}"
            for model_type, model_specs in NS_MODEL_TYPES.items():
                L.info(f"model_type = {model_type}")
                label = f"{planet}_{run_key}_{model_type}"
                ds_dict = {}
                for vrbl in lam_diags.keys():
                    ds_dict[vrbl] = iris.cube.CubeList()
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
                    for vrbl, _calc in lam_diags.items():
                        L.debug(f"vrbl = {vrbl}")
                        input_cl = run.proc
                        if vrbl not in ONLY_LAM:
                            input_cl = input_cl.extract(SS_REGION.constraint)
                        _do_calc(_calc, input_cl, ds_dict[vrbl])

                _dict = {}
                for k, v in ds_dict.items():
                    if len(v) > 0:
                        (cube,) = v.merge()
                        _dict[k] = xr.DataArray.from_iris(cube)
                ds = xr.merge([_dict], compat="override")
                ds.attrs.update(
                    {
                        "region": SS_REGION.to_str(),
                        "model_type": run.model_type,
                        "run_key": run_key,
                        "planet": planet,
                    },
                )
                fname_out = outdir / f"aggr_diags_ns_{run.name}.nc"
                # Write the data to a netCDF file
                ds.to_netcdf(fname_out)
                L.success(f"Saved to {fname_out}")


if __name__ == "__main__":
    t0 = time()
    main()
    L.info(f"Execution time: {time() - t0:.1f}s")
