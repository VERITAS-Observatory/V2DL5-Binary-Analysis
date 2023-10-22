#!/usr/bin/python3
"""
Generate a run list for a given sky region or a given target.

Reads in a yaml file (see example `examples/run_selection.yaml`)
with criteria for selecting runs.

Generates three output files:

1. A list of runs that satisfy the criteria.
2. A detailed printout of the observation table for all runs which satisfy the criteria.
3. A detailed printout of the observation table for all runs which do not satisfy the criteria.

Example:

.. code-block:: console

        $ python generate_runlist.py \
            --obs_table /path/to/ob_table.fits \
            --config examples/run_selection.yaml \
            --output_dir run_lists

"""

import argparse
import logging

import v2dl5.configuration
import v2dl5.run_lists

_logger = logging.getLogger(__name__)


def _parse():
    """
    Parse command line arguments.

    """
    parser = argparse.ArgumentParser(description="Generate a run list.")

    parser.add_argument(
        "--obs_table",
        type=str,
        required=True,
        help="Path to observation table.",
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to configuration file.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory.",
    )

    return parser.parse_args()


def main():
    """
    Generate a run list.

    """
    logging.root.setLevel(logging.INFO)

    args_dict = v2dl5.configuration.configuration(args=_parse(), generate_dqm_run_list=True)

    v2dl5.run_lists.generate_run_list(args_dict)


if __name__ == "__main__":
    main()
