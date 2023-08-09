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
    _plot(plot_name=f"{dataset.name}_{dataset.models[0].name}_fit", output_dir=output_dir)


def plot_flux_points(flux_point_dataset, output_dir=None):
    """
    Plot flux points

    """

    _, ax = plt.subplots()
    flux_point_dataset.plot(ax=ax, sed_type="dnde", color="darkorange")
    flux_point_dataset.plot_ts_profiles(ax=ax, sed_type="dnde")
    _plot(plot_name="flux_points", output_dir=output_dir)


def plot_sed(flux_point_dataset, output_dir):
    """
    Plot spectral energy distribution

    """

    kwargs_model = {"color": "grey", "ls": "--", "sed_type": "dnde"}
    kwargs_fp = {"color": "black", "marker": "o", "sed_type": "dnde"}
    flux_point_dataset.plot_spectrum(kwargs_fp=kwargs_fp, kwargs_model=kwargs_model)
    _plot(plot_name="spectrum", output_dir=output_dir)
    flux_point_dataset.plot_residuals(method="diff/model")
    _plot(plot_name="residuals", output_dir=output_dir)


def plot_light_curve(light_curve, output_dir):
    """
    Plot light curve

    """

    fig, ax = plt.subplots(
        figsize=(8, 6),
        gridspec_kw={"left": 0.16, "bottom": 0.2, "top": 0.98, "right": 0.98},
    )

    light_curve.plot(ax=ax, marker="o", label="per observation")
    _plot(plot_name="light_curve_per_obs", output_dir=output_dir)


def _plot(plot_name=None, output_dir=None):
    """
    Plotting helper function

    """
    if output_dir is not None:
        _ofile = f"{output_dir}/{plot_name}.png"
        logging.info("Plotting %s", _ofile)
        plt.savefig(_ofile)
    else:
        plt.show()
    plt.clf()
