# coding: utf-8
"""Common objects for mean climate diagnostics of global UM runs."""
# Commonly used standard library tools

# Scientific stack
from cf_units import Unit

import iris
from iris.analysis.calculus import _coord_cos, _curl_differentiate, _curl_regrid, differentiate
from iris.analysis.maths import apply_ufunc
from iris.util import reverse
from iris.experimental import stratify

import metpy.calc as metcalc
import metpy.constants as metconst
import metpy.units as metunits

import numpy as np

# My packages and local scripts
from aeolus.calc import (
    last_year_mean,
    integrate,
    ghe_norm,
    heat_redist_eff,
    bond_albedo,
    spatial,
    minmaxdiff,
    meridional_mean,
    region_mean_diff,
    toa_eff_temp,
    toa_cloud_radiative_effect,
    toa_net_energy,
    precip_sum,
    sfc_net_energy,
    sfc_water_balance,
    vertical_cumsum,
    vertical_mean,
    water_path,
    zonal_mean,
)
from aeolus.calc.metpy import preprocess_iris
from aeolus.const import get_planet_radius
from aeolus.const.const import ScalarCube
from aeolus.coord_utils import ensure_bounds, UM_HGT, UM_LATLON, UM_TIME, coord_to_cube
from aeolus.exceptions import MissingCubeError
from aeolus.misc import net_horizontal_flux_to_region
from aeolus.subset import _dim_constr, l_range_constr
from commons import DAYSIDE, NIGHTSIDE


hgt_cnstr_above_5km = l_range_constr(5, 1e99)
hgt_cnstr_0_5km = l_range_constr(0, 5)
hgt_cnstr_0_15km = l_range_constr(0, 15)
hgt_cnstr_5_20km = l_range_constr(5, 20)
hgt_cnstr_0_1km = l_range_constr(0, 1)

MODEL_TIMESTEP = 86400 / 72
FILE_REGEX = r"umglaa.p[a-z]{1}[0]{6}(?P<timestamp>[0-9]{2,4})_00"

ONLY_GLOBAL = ["eta", "b_alb", "t_sfc_diff_dn", "nondim_rossby", "nondim_rhines"]
ONLY_LAM = ["hflux_q", "hflux_t"]

DIAGS = {
    "mse_hdiv": lambda cl: mse_hdiv_mean(cl)["mse"],
    "dse_hdiv": lambda cl: mse_hdiv_mean(cl)["dse"],
    "lse_hdiv": lambda cl: mse_hdiv_mean(cl)["lse"],
    "nondim_rossby": lambda cl: nondim_rossby_deformation_radius(cl),
    "nondim_rhines": lambda cl: nondim_rhines_number(cl.extract(hgt_cnstr_5_20km)),
    "e_net_toa": toa_net_energy,
    "eta": lambda cl: heat_redist_eff(cl, NIGHTSIDE, DAYSIDE),
    "b_alb": bond_albedo,
    "toa_olr": lambda cl: spatial(cl.extract_strict("toa_outgoing_longwave_flux"), "mean"),
    "toa_osr": lambda cl: spatial(cl.extract_strict("toa_outgoing_shortwave_flux"), "mean"),
    "cre_sw": lambda cl: toa_cloud_radiative_effect(cl, kind="sw"),
    "cre_lw": lambda cl: toa_cloud_radiative_effect(cl, kind="lw"),
    "cre_tot": lambda cl: toa_cloud_radiative_effect(cl, kind="total"),
    "gh_norm": ghe_norm,
    "e_net_sfc": sfc_net_energy,
    "t_sfc_diff_dn": lambda cl: region_mean_diff(cl, "surface_temperature", DAYSIDE, NIGHTSIDE),
    "t_sfc": lambda cl: spatial(cl.extract_strict("surface_temperature"), "mean"),
    "t_sfc_min": lambda cl: spatial(cl.extract_strict("surface_temperature"), "min"),
    "t_sfc_max": lambda cl: spatial(cl.extract_strict("surface_temperature"), "max"),
    "t_sfc_extremes": lambda cl: minmaxdiff(cl, name="surface_temperature"),
    "t_eff": toa_eff_temp,
    "wspd_rms": lambda cl: wspd_typical(cl, "rms"),
    "wspd_rms_0_15km": lambda cl: wspd_typical(cl.extract(hgt_cnstr_0_15km), "rms"),
    "wspd_mean": lambda cl: wspd_typical(cl, "mean"),
    "wspd_mean_0_15km": lambda cl: wspd_typical(cl.extract(hgt_cnstr_0_15km), "mean"),
    "wspd_max": lambda cl: wspd_typical(cl, "max"),
    "wspd_max_0_15km": lambda cl: wspd_typical(cl.extract(hgt_cnstr_0_15km), "max"),
    "wvp": lambda cl: spatial(water_path(cl, "water_vapour"), "mean"),
    "lwp": lambda cl: spatial(water_path(cl, "liquid_water"), "mean"),
    "iwp": lambda cl: spatial(water_path(cl, "ice_water"), "mean"),
    "water_balance": sfc_water_balance,
    "precip_total": lambda cl: spatial(precip_sum(cl), "mean"),
    "precip_stra": lambda cl: spatial(precip_sum(cl, ptype="stra"), "mean"),
    "precip_conv": lambda cl: spatial(precip_sum(cl, ptype="conv"), "mean"),
    "precip_rain": lambda cl: spatial(precip_sum(cl, ptype="rain"), "mean"),
    "precip_snow": lambda cl: spatial(precip_sum(cl, ptype="snow"), "mean"),
    "cld_h": lambda cl: spatial(cl.extract_strict("high_type_cloud_area_fraction"), "mean"),
    "cld_m": lambda cl: spatial(cl.extract_strict("medium_type_cloud_area_fraction"), "mean"),
    "cld_l": lambda cl: spatial(cl.extract_strict("low_type_cloud_area_fraction"), "mean"),
    "cloud_frac": lambda cl: spatial(
        cl.extract_strict("cloud_area_fraction_assuming_maximum_random_overlap"), "mean"
    ),
    "vflux_q": lambda cl: mean_vertical_eddy_flux(cl, "specific_humidity"),
    "vflux_t": lambda cl: mean_vertical_eddy_flux(cl, "air_temperature"),
    "t_incr_adv_5_20km": lambda cl: spatial(
        vertical_mean(
            cl.extract_strict(iris.AttributeConstraint(STASH="m01s12i181")).extract(
                hgt_cnstr_5_20km
            ),
            weight_by=cl.extract_strict("air_density").extract(hgt_cnstr_5_20km),
        ),
        "mean",
    ),
    "t_incr_lh_5_20km": lambda cl: spatial(
        vertical_mean(
            latent_heating_rate(cl.extract(hgt_cnstr_5_20km)),
            weight_by=cl.extract_strict("air_density").extract(hgt_cnstr_5_20km),
        ),
        "mean",
    ),
    "hdiv_5_20km": lambda cl: spatial(
        vertical_mean(
            hdiv(cl.extract(hgt_cnstr_5_20km)),
            weight_by=cl.extract_strict("air_density").extract(hgt_cnstr_5_20km),
        ),
        "mean",
    ),
    "hdiv_0_5km": lambda cl: spatial(
        vertical_mean(
            hdiv(cl.extract(hgt_cnstr_0_5km)),
            weight_by=cl.extract_strict("air_density").extract(hgt_cnstr_0_5km),
        ),
        "mean",
    ),
    "dtdz_0_1km": lambda cl: mean_dry_lapse_rate(cl.extract(hgt_cnstr_0_1km)),
    # "hflux_q": lambda cl: total_hflux(
    #     cl.extract(hgt_cnstr_5_20km), "specific_humidity", SS_REGION
    # ),
    # "hflux_t": lambda cl: total_hflux(
    #     cl.extract(hgt_cnstr_5_20km), "air_temperature", SS_REGION
    # ),
}


def air_temperature(cubelist, const=None):
    try:
        cubelist.extract("air_temperature", strict=True)
        return
    except iris.exceptions.ConstraintMismatchError:
        try:
            thta = cubelist.extract("air_potential_temperature", strict=True)
        except iris.exceptions.ConstraintMismatchError:
            raise Exception("")
        if const is None:
            const = thta.attributes["planet_conf"]

        if len(cubelist.extract("dimensionless_exner_function")) == 1:
            exner = cubelist.extract("dimensionless_exner_function", strict=True)
        elif len(cubelist.extract("air_pressure")) == 1:
            pres = cubelist.extract("air_pressure", strict=True)
            exner = (pres / const.reference_surface_pressure.asc) ** (
                const.dry_air_gas_constant / const.dry_air_spec_heat_press
            ).data
        else:
            raise Exception("")
        temp = thta * exner
        temp.rename("air_temperature")
        temp.convert_units("K")
        cubelist.append(temp)


def air_density(cubelist, const=None):
    try:
        cubelist.extract("air_density", strict=True)
        return
    except iris.exceptions.ConstraintMismatchError:
        try:
            temp = cubelist.extract("air_temperature", strict=True)
            pres = cubelist.extract("air_pressure", strict=True)
        except iris.exceptions.ConstraintMismatchError:
            raise Exception("")
        if const is None:
            const = pres.attributes["planet_conf"]
        rho = pres / (const.dry_air_gas_constant.asc * temp)
        rho.rename("air_density")
        rho.convert_units("kg m^-3")
        cubelist.append(rho)


def geopotential_height(cubelist, const=None):
    try:
        g_hgt = cubelist.extract("geopotential_height", strict=True)
        return
    except iris.exceptions.ConstraintMismatchError:
        cube_w_height = cubelist.extract(_dim_constr(UM_HGT, strict=False))[0]
        if const is None:
            const = cube_w_height.attributes["planet_conf"]
        g_hgt = coord_to_cube(cube_w_height, UM_HGT) * const.gravity.asc
        g_hgt.attributes = {k: v for k, v in cube_w_height.attributes.items() if k != "STASH"}
        ensure_bounds(g_hgt, [UM_HGT])
        g_hgt.rename("geopotential_height")
        g_hgt.convert_units("m^2 s^-2")
        cubelist.append(g_hgt)


def calc_derived_cubes(cubelist, const=None):
    """Calculate additional variables."""
    if const is None:
        try:
            const = cubelist[0].attributes["planet_conf"]
        except KeyError:
            raise Exception("")
    air_temperature(cubelist, const=const)
    air_density(cubelist, const=const)
    geopotential_height(cubelist, const=const)


def wspd_typical(cubelist, aggr="median"):
    u = cubelist.extract_strict("x_wind")
    v = cubelist.extract_strict("y_wind")
    return spatial((u ** 2 + v ** 2) ** 0.5, aggr).collapsed(
        [UM_TIME, UM_HGT], getattr(iris.analysis, aggr.upper())
    )


def total_hflux(cubelist, varname, region):
    scalar = cubelist.extract_strict(varname).copy()
    rho = cubelist.extract_strict("air_density")
    scalar *= rho
    scalar.coord(UM_HGT).bounds = None
    scalar.rename(f"{varname}_weighted_by_density")
    u = cubelist.extract_strict("x_wind").copy()
    v = cubelist.extract_strict("y_wind").copy()
    for cube in (scalar, u, v):
        for coord in UM_LATLON:
            cube.coord(coord).bounds = None
    return net_horizontal_flux_to_region(scalar, region, u, v)


def vertical_flux(cubelist, quantity, weight_by_density=True):
    """Vertical flux."""
    w = cubelist.extract_strict("upward_air_velocity")
    q = cubelist.extract_strict(quantity)
    vf = w * q
    if weight_by_density:
        vf *= cubelist.extract_strict("air_density")
    vf.rename(f"vertical_flux_of_{quantity}")
    return vf


def vertical_eddy_flux(cubelist, quantity, w_dens=True):
    """Vertical flux."""
    w = cubelist.extract_strict("upward_air_velocity")
    q = cubelist.extract_strict(quantity)
    w_eddy = w - spatial(w, "mean")
    q_eddy = q - spatial(q, "mean")
    vf = w_eddy * q_eddy
    if w_dens:
        vf *= cubelist.extract_strict("air_density")
    vf.rename(f"vertical_flux_of_{quantity}")
    return vf


def mean_vertical_eddy_flux(cubelist, quantity, weight_by_density=True, coord=UM_HGT):
    vf = spatial(vertical_eddy_flux(cubelist, quantity, w_dens=True), "mean")
    if weight_by_density:
        weight_by = spatial(cubelist.extract_strict("air_density"), "mean")
    else:
        weight_by = None
    mean_vf = vertical_mean(vf, coord=coord, weight_by=weight_by)
    return mean_vf


def meridional_mass_streamfunction(cubelist, z_coord=UM_HGT):
    v = cubelist.extract_strict("y_wind")
    const = v.attributes["planet_conf"]
    if v.coord(z_coord).units.is_convertible("m"):
        rho = cubelist.extract_strict("air_density")
        rho.coord(z_coord).bounds = None
        v.coord(z_coord).bounds = None
        integrand = zonal_mean(v * rho)
        integrand = reverse(integrand, z_coord)
        res = -1 * vertical_cumsum(integrand, coord=z_coord)
        res = reverse(res, z_coord)
    elif v.coord(z_coord).units.is_convertible("Pa"):
        # pres = cubelist.extract_strict("air_pressure")
        # mmstreamf_const /= const.gravity.asc
        raise NotImplementedError()
    coslat = apply_ufunc(np.cos, apply_ufunc(np.deg2rad, coord_to_cube(res, UM_LATLON[0])))
    coslat.units = "1"
    mmstreamf_const = 2 * np.pi * coslat * const.radius.asc
    res = res * mmstreamf_const
    res.rename("meridional_mass_streamfunction")
    res.convert_units("kg s^-1")
    return res


def hdiv(cubelist, i_name_or_constr="x_wind", j_name_or_constr="y_wind"):
    r"""
    Calculate horizontal divergence.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Cubelist containing i-component and j-component of the vector.
    i_name_or_constr: str or iris.Constraint, optional
        Constraint to extract i_cube.
    j_name_or_constr: str or iris.Constraint, optional
        Constraint to extract j_cube.

    Returns
    -------
    div: iris.cube.Cube
        Cube of horizontal divergence.

    Notes
    -----
    Divergence in spherical coordinates is defined as

    .. math::

        \nabla\cdot \vec A = \frac{1}{r cos \theta} (
        \frac{\partial \vec A_\lambda}{\partial \lambda}
        + \frac{\partial}{\partial \theta}
        (\vec A_\theta cos \theta))

    where \lambda is longitude, \theta is latitude.
    """
    i_cube = cubelist.extract_strict(i_name_or_constr)
    j_cube = cubelist.extract_strict(j_name_or_constr)

    x_coord = i_cube.coord(axis="X")
    y_coord = i_cube.coord(axis="Y")

    y_dim = i_cube.coord_dims(y_coord)[0]

    horiz_cs = i_cube.coord_system("CoordSystem")

    # Check for spherical coords
    spherical_coords = isinstance(
        horiz_cs, (iris.coord_systems.GeogCS, iris.coord_systems.RotatedGeogCS)
    )
    if spherical_coords:
        if (y_coord.name() not in ["latitude", "grid_latitude"]) or (
            x_coord.name() not in ["longitude", "grid_longitude"]
        ):
            raise ValueError(
                "Expecting latitude as the y coord and "
                "longitude as the x coord for spherical derivatives."
            )

        # Get the radius of the planet
        r = get_planet_radius(i_cube)
        r_unit = Unit("m")

        lon_coord = x_coord.copy()
        lat_coord = y_coord.copy()
        lon_coord.convert_units("radians")
        lat_coord.convert_units("radians")
        lat_cos_coord = _coord_cos(lat_coord)

        # j-component: \frac{\partial}{\partial \theta} (\vec A_\theta cos \theta))
        temp = iris.analysis.maths.multiply(j_cube, lat_cos_coord, y_dim)
        djcos_dtheta = _curl_differentiate(temp, lat_coord)
        prototype_diff = djcos_dtheta

        # i-component: \frac{\partial \vec A_\lambda}{\partial \lambda}
        d_i_cube_dlambda = _curl_differentiate(i_cube, lon_coord)
        d_i_cube_dlambda = _curl_regrid(d_i_cube_dlambda, prototype_diff)
        new_lat_coord = d_i_cube_dlambda.coord(axis="Y")
        new_lat_cos_coord = _coord_cos(new_lat_coord)
        lat_dim = d_i_cube_dlambda.coord_dims(new_lat_coord)[0]

        # Sum and divide
        div = iris.analysis.maths.divide(
            iris.analysis.maths.add(d_i_cube_dlambda, djcos_dtheta),
            r * new_lat_cos_coord,
            dim=lat_dim,
        )
        div.units /= r_unit
        div = div.regrid(i_cube, iris.analysis.Linear())
    else:
        raise NotImplementedError("Non-spherical coordinates are not implemented yet.")
    return div


def nondim_rossby_deformation_radius(cubelist, const=None, method="direct"):
    r"""
    Estimate the non-dimensional Rossby radius of deformation.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input cubelist.
    const: aeolus.const.const.ConstContainer, optional
        If not given, constants are attempted to be retrieved from
        attributes of a cube in the cube list.
    method: str, optional
        Method of calculation.
        "direct" (default): estimate scale height and BV frequency from air temperature.
        "leconte2013": use isothermal approximation.

    Returns
    -------
    iris.cube.Cube
        Cube with collapsed spatial dimensions.

    References
    ----------
    * Leconte et al. (2013), https://doi.org/10.1051/0004-6361/201321042
      (for `method=leconte2013`)
    .. math::

        \lambda_{Rossby} = \sqrt{\frac{NH}{2\Omega R_p}}
        \approx \sqrt{\frac{R}{c_p^{1/2}} \frac{T^{1/2}}{2\Omega R_p}}

    """
    if const is None:
        const = cubelist[0].attributes["planet_conf"]

    omega = (const.day / (2 * np.pi)) ** (-1)
    double_omega_radius = ScalarCube.from_cube(2 * omega * const.radius)

    if method == "direct":
        rho = cubelist.extract_strict("air_density").copy()
        bv_freq_proxy = spatial(vertical_mean(bv_freq_sq(cubelist), weight_by=rho), "mean") ** 0.5
        temp_proxy = spatial(
            vertical_mean(cubelist.extract_strict("air_temperature"), weight_by=rho), "mean"
        )
        scale_height = const.dry_air_gas_constant.asc * temp_proxy / const.gravity.asc
        nondim_rossby = (bv_freq_proxy * scale_height / double_omega_radius.asc) ** 0.5
    elif method == "leconte2013":
        temp_proxy = toa_eff_temp(cubelist)
        _const_term = ScalarCube.from_cube(const.dry_air_gas_constant / double_omega_radius)
        sqrt_t_over_cp = (temp_proxy / const.dry_air_spec_heat_press.asc) ** 0.5
        nondim_rossby = (sqrt_t_over_cp * _const_term.asc) ** 0.5
    nondim_rossby.convert_units("1")
    nondim_rossby.rename("nondimensional_rossby_deformation_radius")
    return nondim_rossby


def nondim_rhines_number(cubelist, const=None, wspd_aggr="mean"):
    r"""
    Estimate the non-dimensional Rhines number.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Cubelist containing u- and v-components of wind vector.
    const: aeolus.const.const.ConstContainer, optional
        If not given, constants are attempted to be retrieved from
        attributes of a cube in the cube list.

    Returns
    -------
    nondim_rhines: iris.cube.CubeList
        Cube with collapsed spatial dimensions.

    References
    ----------
    * Haqq-Misra et al. (2017), https://doi.org/10.3847/1538-4357/aa9f1f
    .. math::
        \lambda_{Rhines} = \pi \sqrt{U/\beta} / R_p
    """
    # wspd_aggr = "rms"

    def _wspd_typical(cubelist, aggr="mean"):
        u = cubelist.extract_strict("x_wind")
        v = cubelist.extract_strict("y_wind")
        return spatial((u ** 2 + v ** 2) ** 0.5, aggr)

    if const is None:
        const = cubelist[0].attributes["planet_conf"]

    omega = (const.day / (2 * np.pi)) ** (-1)
    double_omega_radius = ScalarCube.from_cube(2 * omega * const.radius)
    u_proxy = _wspd_typical(cubelist, wspd_aggr) / 2
    u_proxy = u_proxy.collapsed(UM_HGT, getattr(iris.analysis, wspd_aggr.upper()))
    nondim_rhines = np.pi * (u_proxy / double_omega_radius.asc) ** 0.5
    nondim_rhines.convert_units("1")
    nondim_rhines.rename("nondimensional_rhines_number")
    return nondim_rhines


def mean_dry_lapse_rate(cubelist, coord=UM_HGT):
    temp = cubelist.extract_strict("air_temperature")
    rho = cubelist.extract_strict("air_density")
    temp = spatial(temp, "mean")
    rho = spatial(rho, "mean")
    dtdz = differentiate(temp, coord)
    dtdz = dtdz.interpolate([(coord, rho.coord(coord).points)], iris.analysis.Linear())
    res = vertical_mean(dtdz, coord=coord, weight_by=rho)
    return res


def moist_static_energy(cubelist, const=None):
    """
    Calculate dry and latent components of the moist static energy.

    .. math::
        h = c_p T + g z + L q
    """
    if const is None:
        const = cubelist[0].attributes["planet_conf"]

    # Get all the necessary cubes
    try:
        temp = cubelist.extract_strict("air_temperature")
        ghgt = cubelist.extract_strict("geopotential_height")
        q = cubelist.extract_strict("specific_humidity")
        q = preprocess_iris(metcalc.mixing_ratio_from_specific_humidity)(q)
    except iris.exceptions.ConstraintMismatchError:
        varnames = ["air_temperature", "geopotential_height", "specific_humidity"]
        raise MissingCubeError(
            f"{varnames} required to calculate mixing ratio are missing from cubelist:\n{cubelist}"
        )
    # Dry component: c_p T + g z
    dry = temp * const.dry_air_spec_heat_press.asc + ghgt
    dry.rename("dry_static_energy")

    # Latent component: L q
    latent = q * const.water_heat_vaporization.asc
    latent.rename("latent_static_energy")

    # Total
    mse = dry + latent
    mse.rename("moist_static_energy")
    return dict(mse=mse, dse=dry, lse=latent)


def mse_hdiv_mmean(cubelist, zcoord=UM_HGT):
    rho = last_year_mean(cubelist.extract_strict("air_density"))
    ensure_bounds(rho, [zcoord])
    u = cubelist.extract_strict("x_wind")
    ensure_bounds(u, [zcoord])
    v = cubelist.extract_strict("y_wind")
    ensure_bounds(v, [zcoord])

    mse_cmpnts = moist_static_energy(cubelist)

    mse_hdiv_cmpnts = {}
    for key, cmpnt in mse_cmpnts.items():
        flux_x = last_year_mean(u * cmpnt)
        flux_x = integrate(rho * flux_x, zcoord)
        flux_x.rename(f"eastward_{cmpnt.name()}_flux")
        flux_y = last_year_mean(v * cmpnt)
        flux_y = integrate(rho * flux_y, zcoord)
        flux_y.rename(f"northward_{cmpnt.name()}_flux")

        fluxes_vec = iris.cube.CubeList([flux_x, flux_y])

        flux_div = hdiv(fluxes_vec, *[i.name() for i in fluxes_vec])

        flux_div_mean = meridional_mean(flux_div)
        flux_div_mean.rename(f"integrated_meridional_mean_flux_divergence_of_{cmpnt.name()}")
        flux_div_mean.convert_units("W m^-2")
        mse_hdiv_cmpnts[key] = flux_div_mean
    return mse_hdiv_cmpnts


def mse_hdiv_mean(cubelist, zcoord=UM_HGT):
    rho = cubelist.extract_strict("air_density")
    ensure_bounds(rho, [zcoord])
    u = cubelist.extract_strict("x_wind")
    ensure_bounds(u, [zcoord])
    v = cubelist.extract_strict("y_wind")
    ensure_bounds(v, [zcoord])

    mse_cmpnts = moist_static_energy(cubelist)

    mse_hdiv_cmpnts = {}
    for key, cmpnt in mse_cmpnts.items():
        ensure_bounds(cmpnt, [zcoord])
        flux_x = u * cmpnt
        flux_x.rename(f"eastward_{cmpnt.name()}_flux")
        flux_y = v * cmpnt
        flux_y.rename(f"northward_{cmpnt.name()}_flux")

        fluxes_vec = iris.cube.CubeList([flux_x, flux_y])

        flux_div = hdiv(fluxes_vec, *[i.name() for i in fluxes_vec])

        flux_div_vmean = integrate(rho * flux_div, zcoord)
        flux_div_mean = spatial(flux_div_vmean, "mean")
        try:
            flux_div_mean = flux_div_mean.collapsed(UM_TIME, iris.analysis.MEAN)
        except iris.exceptions.CoordinateCollapseError:
            pass
        flux_div_mean.rename(f"integrated_horizontal_divergence_of_{cmpnt.name()}")
        flux_div_mean.convert_units("W m^-2")
        mse_hdiv_cmpnts[key] = flux_div_mean
    return mse_hdiv_cmpnts


def latent_heating_rate(cubelist):
    """Temperature increments due to LS precip and (if available) convection."""
    lsppn = cubelist.extract_strict(
        "change_over_time_in_air_temperature_due_to_stratiform_precipitation"
    )
    lh = lsppn.copy()
    try:
        lh += cubelist.extract_strict("change_over_time_in_air_temperature_due_to_convection")
    except iris.exceptions.ConstraintMismatchError:
        pass
    lh.rename("change_over_time_in_air_temperature_due_to_latent_heat_release")
    return lh


def bv_freq_sq(cubelist, const=None):
    """
    Brunt-Vaisala frequency squared (depends on MetPy).

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Cubelist containing potential temperature with a height coordinate.
    const: aeolus.const.const.ConstContainer, optional
        If not given, constants are attempted to be retrieved from
        attributes of a cube in the cube list.

    Returns
    -------
    iris.cube.Cube
        Cube of N^2
    """
    if const is None:
        const = cubelist[0].attributes["planet_conf"]
    theta = cubelist.extract_strict("air_potential_temperature")
    heights = theta.coord(UM_HGT).points * metunits.units.metre
    metconst.g = const.gravity.data * metunits.units(const.gravity.units.format(1))
    res = preprocess_iris(metcalc.brunt_vaisala_frequency_squared)(
        heights, theta, axis=theta.coord_dims(theta.coord(UM_HGT))[0]
    )
    res.rename("brunt_vaisala_frequency_squared")
    res.convert_units("s^-2")
    return res


def interp_to_pressure_levels(
    cubelist, constraint, levels, coord_src=UM_HGT, pressure_cube="air_pressure"
):
    """
    A wrapper for relevel() function.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        List of cubes containing a `varname` cube and a cube of air pressure.
    constraint: str or iris.Constraint
        Variable name or constraint to extract a cube from `cubelist`.
    levels: array-like
        Sequence of pressure levels (same units as `pressure_cube`).

    Returns
    -------
    iris.cube.Cube
        Cube of `varname` interpolated to pressure levels.
    """
    cube = cubelist.extract_strict(constraint)
    pres = cubelist.extract_strict(pressure_cube)
    cube_plev = stratify.relevel(cube, pres, levels, axis=coord_src)
    cube_plev.coord(pressure_cube).attributes = {}
    return iris.util.squeeze(cube_plev)
