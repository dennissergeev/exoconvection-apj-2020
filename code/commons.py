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
# labels += [
#     i.with_suffix("").name.replace("rose-app-", "")
#     for i in sorted(Path(str(topdir) + "_conf").glob("rose-app-*"))
# ]
ENS_LABELS = [
    "base",
    "grcs_ens",
    "amdet_fac_0p01",
    "amdet_fac_0p1",
    "amdet_fac_10p0",
    "amdet_fac_20p0",
    "base",
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
    "rose-app",
]
