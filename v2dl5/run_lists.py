"""
Run list selection from observation table.

"""

import logging

import astropy.table
import astropy.units as u
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

    #    obs_table = _apply_selection_cuts_n_obs(obs_table, n_obs_min=2)
    #    obs_table = _apply_selection_cuts_n_tel(obs_table, n_tel_min=2)
    #    obs_table = _apply_selection_cuts_n_tel_type(obs_table, n_tel_type_min=2)

    return _obs_table


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

    _logger.info("Write run list to %s", output_dir)

    # write a single column of the astropy table to a text file
    column_data = obs_table["OBS_ID"]
    ascii.write(
        column_data,
        f"{output_dir}/run_list.txt",
        overwrite=True,
        format="no_header",
        delimiter="\n",
    )
