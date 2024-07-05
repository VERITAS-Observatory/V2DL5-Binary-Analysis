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
        "--period",
        type=float,
        required=False,
        default=317.3,
        help="Orbital period (in units of days)",
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
    logging.info("orbital period: %.2f" % args.period)
    logging.info("number of bins for averaging: %.1f" % args.orbital_bins)

    binary = {}
    binary["HESS J0632+057"] = {}
    binary["HESS J0632+057"]["name"] = "HESS J0632+057"
    binary["HESS J0632+057"]["orbital_period"] = 317.3
    binary["HESS J0632+057"]["mjd_0"] = 54857.0  # Bongiorno et al 2011

    data_reader = v2dl5.light_curves.data_reader.LightCurveDataReader(
        args.configuration, binary=binary.get("HESS J0632+057")
    )
    data_reader.read_data()
    print("Data reader: ", data_reader.data_dict)

    plotter = v2dl5.light_curves.binary_plotting.BinaryLightCurvePlotter(
        data=data_reader.data_dict,
        config=data_reader.config,
        binary=binary.get("HESS J0632+057"),
    )
    plotter.plot_flux_vs_time("MJD", None, None, args.plot_type)
    plotter.plot_flux_vs_time("orbital phase", None, None, args.plot_type)


#    icrc2019Plots = False
#    donotplotaverage = False
#    lightCurvePlotting.plotLightCurveData(fDataDict,
#                                          PlotInstruments, PlotVariable,
#                                          OrbitalPeriod_BL, OrbitalPeriodBins,
#                                          icrc2019Plots, donotplotaverage)


if __name__ == "__main__":
    main()
