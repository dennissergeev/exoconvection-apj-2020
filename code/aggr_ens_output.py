#!/usr/bin/env python
# coding: utf-8
"""Get mean global climate diagnostics of a multi-config run (ensemble)."""
from pathlib import Path
from time import time
import warnings

# Scientific stack
import iris

import numpy as np

import pandas as pd

# My packages and local scripts
from aeolus.calc import last_year_mean
from aeolus.core import Run
from aeolus.exceptions import MissingCubeError

from commons import (
    DAYSIDE,
    ENS_LABELS,
    NIGHTSIDE,
    PLANET_ALIASES,
    SS_REGION,
)
from gl_diag import (
    DIAGS,
    ONLY_GLOBAL,
    ONLY_LAM,
)
import mypaths
from utils import create_logger


# Global definitions and styles
warnings.filterwarnings("ignore")
L = create_logger(Path(__file__))


@L.catch
def main():
    """Main function."""

    def _do_calc(_calc, cl):
        try:
            cube = _calc(cl)
            try:
                cube = last_year_mean(cube)
            except iris.exceptions.CoordinateCollapseError:
                pass
            value = float(cube.data)
        except (MissingCubeError, iris.exceptions.ConstraintMismatchError) as e:
            L.debug(f"Caught exception:\n{e}")
            value = np.nan
        return value

    for planet in PLANET_ALIASES.keys():
        L.info(f"planet = {planet}")
        topdir = mypaths.sadir / f"{planet}_grcs_ensemble"
        for label in ENS_LABELS:
            L.info(f"label = {label}")
            procdir = topdir / label / "_processed"

            # Initialise a `Run` by loading processed data
            run = Run(
                files=procdir / f"{planet}_{label}.nc",
                name=label,
                planet=planet,
                processed=True,
            )

            data = {}
            for vrbl, _calc in DIAGS.items():
                L.debug(f"vrbl = {vrbl}")
                if vrbl in ONLY_GLOBAL + ONLY_LAM:
                    data[vrbl] = _do_calc(_calc, run.proc)
                else:
                    for _region, suffix in zip(
                        [None, DAYSIDE, NIGHTSIDE, SS_REGION], ["", "_d", "_n", "_ss"]
                    ):
                        input_cl = run.proc
                        if _region is not None:
                            input_cl = input_cl.extract(_region.constraint)
                        data[f"{vrbl}{suffix}"] = _do_calc(_calc, input_cl)

            df = pd.DataFrame(data, index=pd.Index(name="run", data=[run.name])).T
            fname_out = procdir / f"{planet}_grcs_ens_{label}_aggr_diag.csv"
            df.to_csv(fname_out)
            L.success(f"Saved to {fname_out}")


if __name__ == "__main__":
    t0 = time()
    main()
    L.info(f"Execution time: {time() - t0:.1f}s")
