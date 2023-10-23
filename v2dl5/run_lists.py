"""
Run list selection from observation table.

"""

import logging

import astropy.table
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import ascii

_logger = logging.getLogger(__name__)


def generate_run_list(args_dict, target):
    """
    Read observation index table, apply selection cuts and write run list.

    """

    _logger.info("Generate run list. from %s", args_dict["obs_table"])

    obs_table = astropy.table.Table.read(args_dict["obs_table"])

    obs_table = _apply_selection_cuts(obs_table, args_dict, target)

    _logger.info("Selected %d runs.", len(obs_table))

    _write_run_list(obs_table, args_dict["output_dir"])


def _apply_selection_cuts(obs_table, args_dict, target):
    """
    Apply selection cuts to observation table.

    Parameters
    ----------
    obs_table : `~astropy.table.Table`
        Observation table.

    """

    _logger.info("Apply selection cuts.")

    _obs_table = _apply_cut_target(obs_table, args_dict, target)
    _obs_table = _apply_cut_atmosphere(_obs_table, args_dict)
    _obs_table = _apply_cut_dqm(_obs_table, args_dict)
    _obs_table = _apply_cut_ontime_min(_obs_table, args_dict)

    return _obs_table


def _apply_cut_ontime_min(obs_table, args_dict):
    """
    Apply ontime min cut.

    """

    try:
        ontime_min = u.Quantity(args_dict["dqm"]["ontime_min"]).to(u.s)
    except KeyError:
        _logger.error("KeyError: dqm.ontime_min")
        raise

    mask = [row["ONTIME"] > ontime_min.value for row in obs_table]
    _logger.info(f"Remove {mask.count(False)} runs with ontime < {ontime_min}")
    obs_table = obs_table[mask]
    _logger.info(f"Minimum run time: {np.min(obs_table['ONTIME'])} s")

    return obs_table


def _apply_cut_dqm(obs_table, args_dict):
    """
    Apply dqm cuts

    """

    try:
        dqm_stat = args_dict["dqm"]["dqmstat"]
    except KeyError:
        _logger.error("KeyError: dqm.dqmstat")
        raise

    mask = [row["DQMSTAT"] in dqm_stat for row in obs_table]
    _logger.info(f"Remove {mask.count(False)} runs with dqm status not not in {dqm_stat}")
    obs_table = obs_table[mask]
    _logger.info(f"Selected dqm status {np.unique(obs_table['DQMSTAT'])}")

    return obs_table


def _apply_cut_atmosphere(obs_table, args_dict):
    """
    Remove all fields in column "WEATHER" which are not in the list of args_dict.atmosphere.weather

    """

    try:
        weather = args_dict["atmosphere"]["weather"]
    except KeyError:
        _logger.error("KeyError: atmosphere.weather")
        raise

    try:
        mask = [row["WEATHER"][0] in weather for row in obs_table]
    except IndexError:
        _logger.error("IndexError: weather")
        raise
    _logger.info(f"Remove {mask.count(False)} runs with weather not in {weather}")
    obs_table = obs_table[mask]
    _logger.info(f"Selected weather conditions {np.unique(obs_table['WEATHER'])}")

    return obs_table


def _apply_cut_target(obs_table, args_dict, target):
    """
    Apply target cut.

    """

    try:
        obs_cone_radius = u.Quantity(args_dict["observations"]["obs_cone_radius"])
    except KeyError:
        _logger.error("KeyError: observations.obs_cone_radius")
        raise

    angular_separation = target.separation(
        SkyCoord(obs_table["RA_PNT"], obs_table["DEC_PNT"], unit=(u.deg, u.deg))
    )
    obs_table = obs_table[angular_separation < obs_cone_radius]

    _logger.info("Selecting %d runs from observation cone around %s", len(obs_table), target)

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

    _logger.info(f"Write run list to {output_dir}/run_list.txt")

    # write a single column of the astropy table to a text file
    column_data = obs_table["OBS_ID"]
    ascii.write(
        column_data,
        f"{output_dir}/run_list.txt",
        overwrite=True,
        format="no_header",
        delimiter="\n",
    )

    _logger.info(f"Write run table with selected runs to {output_dir}/run_list.fits.gz")

    obs_table.write(f"{output_dir}/run_list.fits.gz", overwrite=True)
