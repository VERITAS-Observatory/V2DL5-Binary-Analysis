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

import v2dl5.target


def _parse():

    parser = argparse.ArgumentParser(
        description="Perform 1D reflected region analysis for a given run list."
    )

    parser.add_argument(
        "--run_list",
        type=str,
        required=True,
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

    logger = logging.getLogger()
    logging.basicConfig(level=logging.INFO)

    args = _parse()

    target = v2dl5.target.Target(
        name=args.target, ra=args.ra, dec=args.dec)

#    runlist = RunList(args.run_list)

#    analysis = Analysis(
#        configuration=args.config, 
#        output_dir=args.output_dir,
#        target=target,
#        runlist=runlist)
#    analysis.run()

if __name__ == "__main__":
    main()
