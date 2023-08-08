"""
Plotting functions for V2DL5
"""

import logging

import matplotlib.pyplot as plt


def plot_fit(dataset, output_dir=None):
    """
    Plot fit results and residuals

    """

    ax_spectrum, _ = dataset.plot_fit()
    ax_spectrum.set_ylim(0.1, 40)
    dataset.plot_masks(ax=ax_spectrum)
    if output_dir is not None:
        _ofile = f"{output_dir}/{dataset.name}_{dataset.models[0].name}_fit.png"
        logging.info("Writing %s fit plot to %s", dataset.name, _ofile)
        plt.savefig(_ofile)
    else:
        plt.show()


def plot_flux_points(flux_points, output_dir=None):
    """
    Plot flux points

    """

    fig, ax = plt.subplots()
    flux_points.plot(ax=ax, sed_type="e2dnde", color="darkorange")
    flux_points.plot_ts_profiles(ax=ax, sed_type="e2dnde")
    if output_dir is not None:
        _ofile = f"{output_dir}/flux_points.png"
        logging.info("Writing flux points fit plot to %s", _ofile)
        plt.savefig(_ofile)
    else:
        plt.show()
