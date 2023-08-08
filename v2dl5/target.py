"""
Target definition.
"""

import logging

from astropy.coordinates import SkyCoord, name_resolve


def get_target(name=None, ra=None, dec=None):
    """
    Defines a SkyCoord object for the target.
    Reads target coordinates from Simbad if target name is given.

    Parameters
    ----------
    name : str
        Target name.
    ra : float
        Target right ascension (deg).
    dec : float
        Target declination (deg).

    """

    _logger = logging.getLogger(__name__)

    target = None
    if name is None:
        target = SkyCoord(ra=ra, dec=dec, unit="deg", frame="icrs")
    else:
        try:
            target = SkyCoord.from_name(name)
            _logger.info("Target %s found in Simbad.", name)
        except name_resolve.NameResolveError:
            _logger.error("Target %s not found in Simbad.", name)
            raise

    _logger.info("Target name: %s", name)
    _logger.info("Target coordinates: %s", target)

    return target
