"""
V2DL5 configuration

Includes reading of configuration from file and default parameters.

"""

import logging
import yaml


def configuration(args):
    """
    V2DL5 configuration.

    Parameters
    ----------
    args : dict
        Command line arguments.

    """

    # Read configuration from file
    args_dict = _default_config()

    args_dict['target'] = args.target
    args_dict['ra'] = args.ra
    args_dict['dec'] = args.dec
    args_dict['run_list'] = args.run_list
    args_dict['output_dir'] = args.output_dir

    return args_dict


def _read_config(config):
    """
    Read configuration from yaml file.

    """

    args_dict = {}
    with open(config, 'r') as stream:
        try:
            args_dict = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    print(args_dict)
    return args_dict


def _default_config():
    """
    Default configuration.

    """

    return {
        'observations': {
            'datastore': '../gammapy',
            'obs_cone': {
                'frame': 'icrs',
                'lon': '83.628700 deg',
                'lat': '22.014700 deg',
                'radius': '5 deg'
            },
            'required_irf': ['aeff', 'edisp']
        },
        'datasets': {
            'type': '1d',
            'stack': False,
            'geom': {
                'axes': {
                    'energy': {
                        'min': '0.05 TeV',
                        'max': '30 TeV',
                        'nbins': 20
                    }, 
                    'energy_true': {
                        'min': '0.05 TeV',
                        'max': '50 TeV',
                        'nbins': 40
                    }
                }
            },
            'on_region': {
                'frame': 'icrs',
                'lon': '83.633 deg',
                'lat': '22.014 deg',
                'radius': '0.08944272 deg'
            },
            'containment_correction': False,
            'safe_mask': {
                'methods': ['aeff-default', 'aeff-max'],
                'parameters': {'aeff_percent': 0.1}},
            'background': {'method': 'reflected'}},
        'fit': {
            'fit_range': {
                'min': '0.1 TeV',
                'max': '20 TeV'
            }, 
            'model': 'pl'
        },
        'flux_points': {
            'energy': {
                'min': '0.1 TeV',
                'max': '20 TeV',
                'nbins': 10
            },
            'source': 'crab'
        }
    }
    