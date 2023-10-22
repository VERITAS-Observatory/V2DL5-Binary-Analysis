"""
Run list selection from observation table.

"""

import logging

import astropy.table
from astropy.io import ascii

_logger = logging.getLogger(__name__)


def generate_run_list(args_dict):
    """
    Read observation index table, apply selection cuts and write run list.

    """

    _logger.info("Generate run list. from %s", args_dict["obs_table"])

    obs_table = astropy.table.Table.read(args_dict["obs_table"])

    obs_table = _apply_selection_cuts(obs_table)

    _write_run_list(obs_table, args_dict["output_dir"])


def _apply_selection_cuts(obs_table):
    """
    Apply selection cuts to observation table.

    Parameters
    ----------
    obs_table : `~astropy.table.Table`
        Observation table.

    """

    _logger.info("Apply selection cuts.")

    #    obs_table = _apply_selection_cuts_n_obs(obs_table, n_obs_min=2)
    #    obs_table = _apply_selection_cuts_n_tel(obs_table, n_tel_min=2)
    #    obs_table = _apply_selection_cuts_n_tel_type(obs_table, n_tel_type_min=2)

    return obs_table


def _write_run_list(obs_table, output_dir):
    """
    Write run list.

    Parameters
    ----------
    obs_table : `~astropy.table.Table`
        Observation table.
    output_dir : str
        Output directory.

    """

    _logger.info("Write run list to %s", output_dir)

    # write a single column of the astropy table to a text file
    column_data = obs_table["OBS_ID"]
    ascii.write(column_data, f"{output_dir}/run_list.txt", overwrite=True, format="no_header")
