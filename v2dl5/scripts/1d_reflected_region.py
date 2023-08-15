#!/usr/bin/python3
"""
Perform 1D reflected region analysis for a given run list.

Example:

.. code-block:: console

    $ python 1D_reflected_region.py \
        --run_list region_run_list.txt \
        --target "Crab" \
        --output_dir 1D_reflected_region_output \
        --config examples/1D_reflected_region_config.yaml

    $ python 1D_reflected_region.py \
        --run_list region_run_list.txt \
        --ra 83.6331 --dec 22.014700 \
        --output_dir 1D_reflected_region_output \
        --config examples/1D_reflected_region_config.yaml

"""

import argparse
import logging

from astropy import units as u

import v2dl5.analysis
import v2dl5.configuration
import v2dl5.data
import v2dl5.sky_regions


def _parse():
    parser = argparse.ArgumentParser(
        description="Perform 1D reflected region analysis for a given run list."
    )

    parser.add_argument(
        "--run_list",
        type=str,
        required=False,
        help="List of run IDs to be processed.",
    )
    parser.add_argument(
        "--target",
        type=str,
        required=False,
        help="Target name.",
    )
    parser.add_argument(
        "--ra",
        type=float,
        required=False,
        help="Target right ascension (deg).",
    )
    parser.add_argument(
        "--dec",
        type=float,
        required=False,
        help="Target declination (deg).",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory.",
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to configuration file.",
    )

    # Get command line arguments
    return parser.parse_args()


def main():
    logging.basicConfig(level=logging.INFO)

    args_dict = v2dl5.configuration.configuration(args=_parse())

    sky_regions = v2dl5.sky_regions.SkyRegions(args_dict=args_dict)

    data = v2dl5.data.Data(
        runlist=args_dict["run_list"],
        data_directory=args_dict["observations"]["datastore"],
        target=sky_regions.target,
        viewcone=args_dict["observations"]["obs_cone"].get("viewcone", 0.5 * u.deg),
    )

    args_dict["datasets"]["on_region"]["radius"] = data.get_on_region_radius()
    sky_regions.get_on_region(args_dict["datasets"]["on_region"])

    analysis = v2dl5.analysis.Analysis(
        args_dict=args_dict,
        sky_regions=sky_regions,
        data=data,
    )
    analysis.run()
    analysis.plot()
    analysis.write()


if __name__ == "__main__":
    main()
