import logging

from astropy.coordinates import SkyCoord, Angle, name_resolve
from astroquery.simbad import Simbad


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
        except name_resolve.NameResolveError:
            _logger.error(f"Target {name} not found in Simbad.")
            raise

    _logger.info(f"Target name: {name}")
    _logger.info(f"Target coordinates: {target}")

    return target
