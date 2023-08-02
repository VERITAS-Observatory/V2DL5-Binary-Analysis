from astropy.coordinates import SkyCoord, Angle
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

        self.name = name
        self.ra = ra
        self.dec = dec

        if self.name is not None:
            self._get_coordinates()

    def _get_coordinates(self):
        """
        Get target coordinates from Simbad.
        """


        result_table = Simbad.query_object(self.name)

        if result_table is None:
            raise ValueError("Target name not found in Simbad.")

        self.ra = result_table["RA_d"][0]
        self.dec = result_table["DEC_d"][0]