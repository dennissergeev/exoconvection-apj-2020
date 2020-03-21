"""Common plotting functions and parameters."""
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point

from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

import numpy as np

from scipy import interpolate

from aeolus.coord_utils import UM_LATLON
from aeolus.plot import GeoAxesGrid, label_global_map_gridlines
from aeolus.util import subplot_label_generator


CART_KW = dict(transform=ccrs.PlateCarree())

MARKER_KW = dict(
    trap1e=dict(marker="X"),
    proxb=dict(marker="o"),
    grcs=dict(color="C0"),
    llcs_all_rain=dict(color="C1"),
    acoff=dict(color="C2"),
)

# Locations of grid lines on maps
XLOCS = np.arange(-180, 181, 90)
YLOCS = np.arange(-90, 91, 30)


def add_aux_yticks(
    ax, src_points, src_values, target_points, twin_ax_ylim=None, twin_ax_inv=False
):
    """
    Add Y-axis ticks at desired locations.

    Examples
    --------
    >>> ax = plt.axes()
    >>> ax.plot(temperature, pressure)
    >>> add_aux_yticks(ax, heights, pressure, [0, 10, 40], twin_ax_ylim=[0, 40], twin_ax_inv=True)
    """
    int_func = interpolate.interp1d(src_points, src_values, fill_value="extrapolate")
    new_points = int_func(target_points)
    _ax = ax.twinx()
    _ax.set_ylim(twin_ax_ylim)
    _ax.set_yticks(new_points)
    _ax.set_yticklabels(target_points)
    if twin_ax_inv:
        _ax.invert_yaxis()
    return _ax


def add_custom_legend(ax, styles_and_labels, **leg_kw):
    """
    Add a custom legend to matplotlib axes.

    Parameters
    ----------
    ax: matplotlib.axes._subplots.AxesSubplot
        Axes where to put the legend.
    styles_and_labels: dict
        Dictionary with labels as keys and a dictionary of plot
        keywords as values.
    leg_kw: dict, optional
        Keyword arguments passed to `legend()` function.

    Example
    -------
    >>> import matplotlib.pyplot as plt
    >>> ax = plt.axes()
    >>> my_dict = dict(foo=dict(color='C0', marker="X"),
                       bar=dict(color='C1', marker="o"))
    >>> add_custom_legend(ax, my_dict, loc=2, title="blah")

    """
    lines = [Line2D([0], [0], **style) for style in styles_and_labels.values()]
    leg = ax.legend(lines, styles_and_labels.keys(), **leg_kw)
    if ax.legend_ is not None:
        ax.add_artist(leg)


def use_style():
    """Load custom matplotlib style sheet."""
    plt.style.use("paper.mplstyle")


def make_map_figure(ncols, nrows, **axgr_kw):
    """
    Make a figure with a grid of cartopy axes with the Robinson projection.

    Parameters
    ----------
    ncols: int
        Number of columns
    nrows: int
        Number of rows
    axgr_kw: dict, optional
        Parameters passed to `aeolus.plot.cart.GeoAxesGrid`.

    Returns
    -------
    matplotlib.figure.Figure, aeolus.plot.cart.GeoAxesGrid
        The figure and axes grid.
    """

    iletters = subplot_label_generator()

    fig = plt.figure(figsize=(8 * ncols, 4 * nrows))

    axgr = GeoAxesGrid(
        fig, 111, projection=ccrs.Robinson(), nrows_ncols=(nrows, ncols), **axgr_kw
    )
    for ax in axgr.axes_all:
        label_global_map_gridlines(
            fig, ax, XLOCS[1:-1], YLOCS[1:-1], degree=True, size="x-small", xoff=-15
        )
        ax.gridlines(xlocs=XLOCS, ylocs=YLOCS, crs=ccrs.PlateCarree())
        ax.set_title(f"({next(iletters)})", fontsize="small", pad=5, loc="left")

    return fig, axgr


def draw_scalar_cube(
    cube,
    ax,
    method="contourf",
    cax=None,
    tex_units=None,
    cbar_ticks=None,
    use_cyclic=True,
    **plt_kw,
):
    """
    Plot a cube on a map.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube to plot.
    ax: matplotlib.axes._subplots.AxesSubplot
        Cartopy axes.
    method: str, optional
        Method of plotting, e.g. "contour", "pcolormesh", etc.
    cax: matplotlib.axes._subplots.AxesSubplot or similar
        Axes for the colorbar.
    tex_units: str, optional
        TeX string of cube units to be attached to the colorbar.
    cbar_ticks: sequence, optional
        Colorbar ticks.
    use_cyclic: bool, optional
        Use `cartopy.utils.add_cyclic_point` for the data.
    plt_kw: dict, optional
        Keywords for the plotting method.

    Returns
    -------
    Result of the plotting method.
    """
    lons = cube.coord(UM_LATLON[1]).points
    lats = cube.coord(UM_LATLON[0]).points
    cb_ttl_kw = dict(fontsize="x-small", pad=5)

    if use_cyclic:
        _data, _lons = add_cyclic_point(cube.data, coord=lons)
    else:
        _data, _lons = cube.data, lons
    h = getattr(ax, method)(_lons, lats, _data, **plt_kw, **CART_KW)
    if cax is not None:
        cb = ax.figure.colorbar(h, cax=cax)
        cb_ttl_kw = dict(fontsize="x-small", pad=5)
        if tex_units is not None:
            cb.ax.set_title(f"[{tex_units}]", **cb_ttl_kw)
        if cbar_ticks is not None:
            cb.set_ticks(cbar_ticks)
    return h


def draw_vector_cubes(
    u,
    v,
    ax,
    cax=None,
    tex_units=None,
    cbar_ticks=None,
    mag=(),
    xstride=1,
    ystride=1,
    add_wind_contours=False,
    qk_ref_wspd=None,
    **plt_kw,
):
    """
    Plot vectors of two cubes on a map.

    Parameters
    ----------
    u: iris.cube.Cube
        X-component of the vector.
    v: iris.cube.Cube
        Y-component of the vector.
    ax: matplotlib.axes._subplots.AxesSubplot
        Cartopy axes.
    cax: matplotlib.axes._subplots.AxesSubplot or similar
        Axes for the colorbar.
    tex_units: str, optional
        TeX string of cube units to be attached to the colorbar.
    cbar_ticks: sequence, optional
        Colorbar ticks.
    mag: tuple, optional
        Tuple of numpy arrays to be used for colour-coding the vectors.
    xstride: int, optional
        Stride x-component data.
    ystride: int, optional
        Stride y-component data.
    add_wind_contours: bool, optional
        Add contours of the vector magnitude (wind speed).
    qk_ref_wspd: float, optional
        Reference vector magnitude (wind speed).
        If given, a reference arrow (quiver key) is added to the figure.
    plt_kw: dict, optional
        Keywords passed to quiver().
    """
    cb_ttl_kw = dict(fontsize="x-small", pad=5)
    xsl = slice(xstride, -xstride, xstride)
    ysl = slice(ystride, -ystride, ystride)

    lons = u.coord(UM_LATLON[1]).points
    lats = u.coord(UM_LATLON[0]).points

    h = ax.quiver(
        lons[xsl], lats[ysl], u.data[ysl, xsl], v.data[ysl, xsl], *mag, **plt_kw
    )
    if cax is not None and mag:
        cb = ax.figure.colorbar(h, cax=cax)
        if tex_units is not None:
            cb.ax.set_title(f"[{tex_units}]", **cb_ttl_kw)
        if cbar_ticks is not None:
            cb.set_ticks(cbar_ticks)

    if qk_ref_wspd is not None:
        ax.quiverkey(
            h,
            0.17,
            0.87,
            qk_ref_wspd,
            fr"${qk_ref_wspd}$" + r" $m$ $s^{-1}$",
            labelpos="N",
            labelsep=0.05,
            coordinates="figure",
            color="#444444",
            fontproperties=dict(size="small"),
        )

    if add_wind_contours:
        wspd = (u ** 2 + v ** 2) ** 0.5
        ax.contour(
            lons,
            lats,
            wspd.data,
            transform=plt_kw["transform"],
            levels=np.arange(30, 105, 5),
            cmap="Greens",
        )
