"""Definitions and objects commonly used between scripts."""
from datetime import datetime
from pathlib import Path

import iris

from aeolus.region import Region

# Selected variables
SINGLE_LEVEL_VARS = [
    "surface_temperature",
    "toa_incoming_shortwave_flux",
    "toa_outgoing_longwave_flux",
    "toa_outgoing_longwave_flux_assuming_clear_sky",
    "toa_outgoing_shortwave_flux",
    "toa_outgoing_shortwave_flux_assuming_clear_sky",
    "surface_downwelling_shortwave_flux_in_air",
    "upwelling_shortwave_flux_in_air",
    "surface_downwelling_longwave_flux_in_air",
    "upwelling_longwave_flux_in_air",
    "surface_upward_sensible_heat_flux",
    "surface_upward_latent_heat_flux",
    "convective_rainfall_flux",
    "convective_snowfall_flux",
    "high_type_cloud_area_fraction",
    "low_type_cloud_area_fraction",
    "medium_type_cloud_area_fraction",
    "stratiform_rainfall_flux",
    "stratiform_snowfall_flux",
]
MULTI_LEVEL_VARS = [
    "air_potential_temperature",
    "air_pressure",
    "mass_fraction_of_cloud_liquid_water_in_air",
    "mass_fraction_of_cloud_ice_in_air",
    "relative_humidity",
    "specific_humidity",
    "upward_air_velocity",
]
HORIZ_WINDS_NAMES = ["x_wind", "y_wind"]
HORIZ_WINDS_STASH = ["m01s30i001", "m01s30i002"]
INCR_CONSTR = iris.AttributeConstraint(STASH=lambda x: x.item in [181, 182, 233])

# Common parameters
DAYSIDE = Region(-90, 90, -90, 90, "dayside")
NIGHTSIDE = Region(90, -90, -90, 90, "nightside")

DT_FMT = "%Y%m%dT%H%MZ"

RUN_ALIASES = {"grcs": "MassFlux", "llcs_all_rain": "Adjust", "acoff": "NoCnvPm"}
PLANET_ALIASES = {"trap1e": "Trappist-1e", "proxb": "Proxima b"}
OUTPUT_NAME_PREFIX = (
    f"{'_'.join(PLANET_ALIASES.keys())}__{'_'.join(RUN_ALIASES.keys())}"
)

# Global simulation parameters
GLM_MODEL_TIMESTEP = 86400 / 72
GLM_RUNID = r"umglaa"  # file prefix
GLM_FILE_REGEX = GLM_RUNID + r".p[b,c,d,e]{1}[0]{6}(?P<timestamp>[0-9]{2,4})_00"
GLM_START_DAY = 1440  # by default, use only the last 360 days of a 5 year run

# Nesting Suite parameters
# NS domain
ns_centre = (10, 0)
SS_REGION = Region(
    ns_centre[0] - 10,
    ns_centre[0] + 10,
    ns_centre[1] - 10,
    ns_centre[1] + 10,
    "substellar_region",
)
NS_RUNID = r"umnsaa"
# Simulation time (the UM simulation length is defined using real Earth dates)
NS_CYCLE_TIMES = {
    "trap1e_grcs": {"start": datetime(2009, 7, 27, 9, 0), "ndays": 10},
    "proxb_grcs": {"start": datetime(2009, 7, 27, 9, 0), "ndays": 10},
}
NS_RUN_ALIASES = {"grcs": "MassFlux"}
NS_OUTPUT_NAME_PREFIX = (
    f"{'_'.join(PLANET_ALIASES.keys())}__{'_'.join(NS_RUN_ALIASES.keys())}"
)
NS_MODEL_TYPES = {
    "global": {
        "path": Path("glm") / "um" / f"{GLM_RUNID}_p[a,b,c,d]*",
        "timestep": GLM_MODEL_TIMESTEP,
    },
    "lam": {
        "path": Path("regn_0N10E") / "resn_1" / "ra1t" / "um" / f"{NS_RUNID}_p[b,c,d]*",
        "timestep": 150,
    },
}
NS_COLORS = {
    "grcs": {"global": "deepskyblue", "lam": "navy"},
}

# Ensemble simulation labels
# labels += [
#     i.with_suffix("").name.replace("rose-app-", "")
#     for i in sorted(Path(str(topdir) + "_conf").glob("rose-app-*"))
# ]
ENS_LABELS = [
    "base",
    "amdet_fac_0p01",
    "amdet_fac_0p1",
    "amdet_fac_10p0",
    "amdet_fac_20p0",
    "cape_timescale_14400",
    "cape_timescale_1800",
    "cape_timescale_28800",
    "cape_timescale_300",
    "cnv_cold_pools_1",
    "cnv_cold_pools_2",
    "cnv_wat_load_opt_1",
    "ent_fac_dp_0p01",
    "ent_fac_dp_0p1",
    "ent_fac_dp_0p5",
    "ent_fac_dp_2p0",
    "ent_fac_md_0p01",
    "ent_fac_md_0p1",
    "ent_fac_md_0p5",
    "ent_fac_md_2p0",
    "ent_opt_dp_0",
    "ent_opt_dp_1",
    "ent_opt_dp_2",
    "ent_opt_dp_4",
    "ent_opt_dp_5",
    "ent_opt_dp_6",
    "ent_opt_dp_7",
    "fac_qsat_0p01",
    "fac_qsat_0p1",
    "fac_qsat_5p0",
    "l_mom_0",
    "l_mom_dd_1",
    "mid_cnv_pmin_1000p0",
    "mparwtr_0p0001",
    "mparwtr_0p01",
    "mparwtr_1e-05",
    "qlmin_0p001",
    "qlmin_1e-06",
    "qlmin_1e-08",
    "r_det_0p1",
    "r_det_0p3",
    "r_det_0p6",
    "r_det_1p0",
    "w_cape_limit_0p01",
    "w_cape_limit_0p1",
    "w_cape_limit_10p0",
    "w_cape_limit_1p0",
]
