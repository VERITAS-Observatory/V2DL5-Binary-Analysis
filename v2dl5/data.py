import logging

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from gammapy.data import DataStore


class Data():
    """
    Data class holding data store and observations

    Allows to select data from run list or based on
    target coordinates (and viewcone).

    Parameters
    ----------
    runlist : str
        Path to run list.
    data_directory : str
        Path to data directory.
    ra : float
        Target right ascension (deg).
    dec : float
        Target declination (deg).
    viewcone : float
        Viewcone radius (deg).

    """

    def __init__(
            self,
            runlist = None,
            data_directory = "../gammapy",
            ra = None,
            dec = None,
            viewcone=0.5):
        """
        Initialize Data object.

        Uses runlist if not set to None, otherwise selects data
        according to target coordinates and viewcone.

        """

        self._logger = logging.getLogger(__name__)

        self._data_store = DataStore.from_dir(data_directory)
        if runlist is None:
            self.runs = self._from_target(ra, dec, viewcone)
        else:
            self.runs = self._from_runlist(runlist)

    def get_data_store(self):
        """
        Return data store.

        """

        return self._data_store

    def get_observations(self):
        """
        Return observations.

        Reflected region analysis requires effective area and energy dispersion only.

        """

        available_irf = ["aeff", "edisp"]

        return self._data_store.get_observations(self.runs, required_irf=available_irf)

    def _from_runlist(self, runlist):
        """
        Read runlist from file and select data
        
        Parameters
        ----------
        runlist : str
            Path to run list.

        """

        if runlist is None:
            return None

        try:
            _runs = np.loadtxt(runlist, dtype=int)
        except OSError:
            self._logger.error(f"Run list {runlist} not found.")
            raise

        self._logger.info(f"Reading run list with {len(_runs)} runs from {runlist}")
        return _runs

    def _from_target(self, ra, dec, viewcone):
        """
        Select data based on target coordinates and viewcone.
        
        Parameters
        ----------
        ra : float
            Target right ascension (deg).
        dec : float
            Target declination (deg).
        viewcone : float
            Viewcone radius (deg).

        """

        observations = self._data_store.obs_table
        target = SkyCoord(ra=ra, dec=dec, unit="deg", frame="icrs")
        mask = target.separation(observations.pointing_radec) < viewcone * u.deg
        _runs = observations[mask]["OBS_ID"].data

        self._logger.info(f"Selecting {len(_runs)} runs from viewcone around {target}")
        self._logger.info("WARNING - this is not tested")
        return _runs

