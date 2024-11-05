"""Utilities for plotting."""

import logging
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)


def get_plotting_variable(plot_variable, index):
    """Return string with plotting variable."""
    # plotting variable
    i_plotValue = "flux"
    i_plotError = "flux_err"
    if plot_variable and index < len(plot_variable):
        i_plotValue = plot_variable[index]
        i_plotError = plot_variable[index] + "_err"

    return i_plotValue, i_plotError


def align_yaxis(ax1, v1, ax2, v2):
    """
    Adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1.

    from https://stackoverflow.com/questions/10481990/matplotlib-axis-with-two-scales-shared-origin
    """
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1 - y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny + dy, maxy + dy)


def plot_uncertainty_band(ax, per_stdx):
    """Plot a vertical uncertainty band."""
    ax.axvspan(per_stdx[0], per_stdx[1], alpha=0.5, color="g")
    # ax.axvspan(per_stdx[0], per_stdx[1], alpha=0.5, color='b')


def get_line_width():
    """Return default line width."""
    return 0.5


def get_marker_size():
    """Return default marker size."""
    return 4


def get_color_list(n_colors, plot_type=None):
    """
    Return a list of colors used in all plots.

    Parameters
    ----------
    n_colors: int
        Number of colors
    plot_type: str
        Type of plot

    Returns
    -------
    list
        List of colors

    """
    if plot_type == "icrc2019Plots":
        colormap = plt.cm.Dark2
    elif plot_type == "xray":
        colormap = plt.cm.tab10
    else:
        colormap = plt.cm.tab20b
    return [colormap(i) for i in np.linspace(0, 1, n_colors + 1)]


def get_marker_list():
    """Return a list of markers."""
    return ["o", "s", "P", "^", "D"]


def print_figure(print_name, file_type=".pdf", figure_dir="./figures/"):
    """
    Print a figure into the given format.

    Parameters
    ----------
    print_name :  str
        Name of the figure.
    file_type :   str
        File type of the figure.
    figure_dir :  str
        Directory where the figure is saved.

    """
    # don't allow spaces in print string
    print_name = print_name.replace(" ", "-")
    figure_path = Path(figure_dir)
    figure_path.mkdir(parents=True, exist_ok=True)

    figure_file = figure_path / (print_name + file_type)
    logger.info(f"Printing figure to {figure_file}")
    if "png" in file_type:
        plt.savefig(figure_file, dpi=(400))
    else:
        plt.savefig(figure_file)


def get_orbital_phase_axis_string(orbital_period):
    """Return a string for the orbital phase axis."""
    # strOr = "orbital phase (%.1f d)" % orbital_period
    return "orbital phase"


def get_flux_axis_string(instrument, plot_variable=None, scale_factor=None, energy_flux=False):
    """Flux axis string with the correct unit depending on the type of instrument."""
    if not scale_factor:
        scale_factor = ""

    if "VERITAS" in instrument or "HESS" in instrument or "MAGIC" in instrument:
        if energy_flux:
            pAxisFluxSTR = "Flux E>350 GeV (" + scale_factor + "erg/(cm$^2$ s))"
        else:
            pAxisFluxSTR = "Flux E>350 GeV (" + scale_factor + "1/(cm$^2$s))"
    elif instrument == "XRT":
        pAxisFluxSTR = "Swift rate (0.3-10 keV) (counts/s)"
    elif instrument.find("Optical") >= 0:
        if plot_variable and len(plot_variable) == 1:
            if plot_variable[0] == "ew":
                pAxisFluxSTR = "EW ($\\AA$)"
            elif plot_variable[0] == "vc":
                pAxisFluxSTR = "v$_c$ (km/s)"
            elif plot_variable[0] == "fwhm":
                pAxisFluxSTR = "FWHM (km/s)"
    else:
        pAxisFluxSTR = "Flux 0.3-10 keV (" + scale_factor + "erg/(cm$^2$s))"

    return pAxisFluxSTR


def get_paper_figures_parameters(width=None, height=None, xtick_top=True, legend_font_size=8):
    """Return default paper figure parameters."""
    if width is None:
        width = 6

    if height is None:
        golden_mean = (math.sqrt(5) - 1.0) / 2.0
        height = width * golden_mean

    return {
        "axes.labelsize": 10,  # font size for x and y labels (was 10)
        "axes.titlesize": 8,
        "font.size": 8,  # was 10
        "legend.fontsize": legend_font_size,  # was 10
        "xtick.labelsize": 8,
        "xtick.direction": "in",
        "ytick.labelsize": 8,
        "ytick.direction": "in",
        "xtick.top": xtick_top,
        "ytick.right": True,
        "figure.figsize": [width, height],
        "figure.autolayout": True,
        "font.family": "serif",
    }


def paper_figures(width=None, height=None, columns=1, xtick_top=True):
    """Set figures parameters for paper plotting."""
    paper_parameters = get_paper_figures_parameters(width, height, xtick_top)
    mpl.rcParams.update(paper_parameters)
    width = paper_parameters["figure.figsize"][0] if width is None else width
    height = paper_parameters["figure.figsize"][1] if height is None else height

    plt.gcf().clear()
    plt.figure(figsize=(width, height))

    return plt.gca()
