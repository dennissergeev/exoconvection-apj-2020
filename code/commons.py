"""Definitions and objects commonly used between scripts."""
from aeolus.region import Region


RUN_ALIASES = {"grcs": "MassFlux", "llcs_all_rain": "Adjust", "acoff": "NoCnvPm"}
PLANET_ALIASES = {"trap1e": "Trappist-1e", "proxb": "Proxima b"}
output_name_prefix = (
    f"{'_'.join(PLANET_ALIASES.keys())}__{'_'.join(RUN_ALIASES.keys())}"
)

# Global model parameters
GLM_MODEL_TIMESTEP = 86400 / 72
GLM_RUNID = r"umglaa"  # file prefix
GLM_FILE_REGEX = GLM_RUNID + r".p[b,c,d,e]{1}[0]{6}(?P<timestamp>[0-9]{2,4})_00"
GLM_START_DAY = 1440  # by default, use only the last 360 days of a 5 year run

model_colors = {
    "grcs": {"global": "deepskyblue", "lam": "navy"},
}

DT_FMT = "%Y%m%dT%H%MZ"

DAYSIDE = Region(-90, 90, -90, 90, "dayside")
NIGHTSIDE = Region(90, -90, -90, 90, "nightside")

# NS domain
ns_centre = (10, 0)
SS_REGION = Region(
    ns_centre[0] - 10,
    ns_centre[0] + 10,
    ns_centre[1] - 10,
    ns_centre[1] + 10,
    "substellar_region",
)
