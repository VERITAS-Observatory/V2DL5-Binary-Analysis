"""
V2DL5 configuration

Read configuration from file (or us default parameters).

"""

import logging

import yaml

_logger = logging.getLogger(__name__)


def configuration(args):
    """
    V2DL5 configuration.

    Parameters
    ----------
    args : dict
        Command line arguments.

    """

    args_dict = _default_config()
    if args.config is not None:
        args_dict.update(_read_config_from_file(args.config))

    args_dict["run_list"] = args.run_list
    args_dict["output_dir"] = args.output_dir

    return args_dict


def _read_config_from_file(config):
    """
    Read configuration from yaml file.

    Parameters
    ----------
    config : str
        Path to configuration file.

    """

    _logger.info("Reading configuration from %s", config)

    args_dict = {}
    with open(config, "r", encoding="utf-8") as stream:
        try:
            args_dict = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            _logger.error(exc)
            raise

    _logger.info(args_dict)
    return args_dict


def _default_config():
    """
    Default configuration.

    """

    return {
        "observations": {
            "datastore": "../../../VTS/DL3/v490/point-like/",
            "obs_cone_radius": "5. deg",
            "required_irf": ["aeff", "edisp"],
        },
        "datasets": {
            "type": "1d",
            "stack": False,
            "geom": {
                "axes": {
                    "energy": {"min": "0.05 TeV", "max": "30 TeV", "nbins": 20},
                    "energy_true": {"min": "0.05 TeV", "max": "50 TeV", "nbins": 40},
                }
            },
            "on_region": {
                "target": None,
                "frame": "icrs",
                "lon": "83.633 deg",
                "lat": "22.014 deg",
                "radius": "0.08944272 deg",
            },
            "exclusion_region": {
                "on_radius": "0.5 deg",
                "magnitude_B": 7.0,
                "star_exclusion_radius": "0.3 deg",
                "fov": "5 deg",
                "star_file": "./data/hip_mag9.fits.gz",
            },
            "containment_correction": False,
            "safe_mask": {
                "methods": ["aeff-default", "aeff-max"],
                "parameters": {"aeff_percent": 0.1},
            },
            "background": {"method": "reflected"},
        },
        "fit": {"fit_range": {"min": "0.1 TeV", "max": "20 TeV"}, "model": "pl"},
        "flux_points": {"energy": {"min": "0.1 TeV", "max": "20 TeV", "nbins": 10}, "source": None},
        "light_curve": {
            "energy": {"min": "0.1 TeV", "max": "100 TeV"},
            "time_zone": -7,
        },
    }
