"""Common plotting functions and parameters."""
import cartopy.crs as ccrs

from matplotlib.lines import Line2D
import matplotlib.pyplot as plt


CART_KW = dict(transform=ccrs.PlateCarree())

MARKER_KW = dict(
    trap1e=dict(marker="X"),
    proxb=dict(marker="o"),
    grcs=dict(color="C0"),
    llcs_all_rain=dict(color="C1"),
    acoff=dict(color="C2"),
)


def use_style():
    """Load custom matplotlib style sheet."""
    plt.style.use("paper.mplstyle")


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
