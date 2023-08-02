import logging

from astropy.coordinates import SkyCoord, Angle, name_resolve
from astroquery.simbad import Simbad


class Target():
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

    def __init__(self, name=None, ra=None, dec=None):
        """
        Initialize Target object.
        """

        self._logger = logging.getLogger(__name__)

        self.name = name
        self.target = None
        if name is None:
            self.target = SkyCoord(ra=ra, dec=dec, unit="deg", frame="icrs")
        else:
            try:
                self.target = SkyCoord.from_name(name)
            except name_resolve.NameResolveError:
                self._logger.error(f"Target {name} not found in Simbad.")
                raise

        self._logger.info(f"Target name: {self.name}")
        self._logger.info(f"Target coordinates: {self.target}")
