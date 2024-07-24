#!/usr/bin/python
""""
Plotting of binary light curves.

Includes the following plots

- Light curve vs time
- Light curve vs orbital phase
- Light curve vs orbital phase for each orbit

"""

import argparse
import logging

import v2dl5.binaries as binaries
import v2dl5.light_curves.binary_plotting
import v2dl5.light_curves.data_reader


def _parse():
    """
    Parse command line arguments.

    Returns
    -------
    dict
        Command line parameters.

    """
    parser = argparse.ArgumentParser(description="Binary light-curve plotting.")

    parser.add_argument(
        "--instrument",
        type=str,
        required=True,
        help="Source of data (VERITAS, HESS, Gamma, Swift XRT).",
    )
    parser.add_argument(
        "--configuration",
        type=str,
        required=True,
        help="Configuration describing data files and plotting",
    )
    parser.add_argument(
        "--binary_name",
        type=str,
        required=True,
        default="LS I +61 303",
        help="Binary name (e.g., LS I +61 303; see v2dl5.binaries for definition).",
    )
    parser.add_argument(
        "--orbital_bins",
        type=float,
        required=False,
        default=0.0,
        help="Number of bins in orbital period for averaging.",
    )
    parser.add_argument(
        "--plot_type",
        type=str,
        required=False,
        default=".pdf",
        help="File type for plots (e.g., '.pdf', '.png')",
    )
    parser.add_argument(
        "--figure_dir",
        type=str,
        required=False,
        default="./figures/",
        help="Directory for saving figures",
    )

    return parser.parse_args()


def main():
    """Binary light-curve plotting."""
    logging.root.setLevel(logging.INFO)

    args = _parse()
    if args.orbital_bins < 1:
        if args.instrument in ["VERITAS", "Gamma", "HESS"]:
            args.orbital_bins = 10.0
        else:
            args.orbital_bins = 20.0

    logging.info("Light Curve Analysis - run parameters")
    logging.info("instrument: %s" % args.instrument)
    logging.info("instrument list: %s" % args.configuration)
    logging.info("number of bins for averaging: %.1f" % args.orbital_bins)

    try:
        binary = binaries.binary_properties()[args.binary_name]
    except KeyError:
        raise KeyError(f"Binary {args.binary_name} not found in binaries.py")

    data_reader = v2dl5.light_curves.data_reader.LightCurveDataReader(
        args.configuration, binary=binary
    )
    data_reader.read_data()

    plotter = v2dl5.light_curves.binary_plotting.BinaryLightCurvePlotter(
        data=data_reader.data_dict, config=data_reader.config, binary=binary
    )
    plotter.plot_flux_vs_time(
        "MJD", None, None, file_type=args.plot_type, figure_dir=args.figure_dir
    )
    plotter.plot_flux_vs_time(
        "orbital phase", None, None, file_type=args.plot_type, figure_dir=args.figure_dir
    )
    plotter.plot_flux_vs_phase_for_individual_orbits(
        "VERITAS (anasum)", file_type=args.plot_type, figure_dir=args.figure_dir
    )


if __name__ == "__main__":
    main()
