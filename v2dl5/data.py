"""
Data class holding data store and observations.
"""

import logging

import numpy as np
from astropy import units as u
from gammapy.data import DataStore


class Data:
    """
    Data class holding data store and observations

    Allows to select data from run list or based on
    target coordinates (and viewcone).

    Parameters
    ----------
    runlist : str
        Path to run list.
    data_directory : str
        Path to data directory (holding hdu-index.fits.gz and obs-index.fits.gz).
    target : SkyCoord
        Target coordinates.
    viewcone : float
        Viewcone radius (deg).

    """

    def __init__(self, runlist=None, data_directory=None, target=None, viewcone=0.5 * u.deg):
        """
        Initialize Data object.

        Uses runlist if not set to None, otherwise selects data
        according to target coordinates and viewcone.

        """

        self._logger = logging.getLogger(__name__)

        self._logger.info("Initializing data object from %s", data_directory)
        self._data_store = DataStore.from_dir(data_directory)
        if runlist is None:
            self.runs = self._from_target(target, viewcone)
        else:
            self.runs = self._from_runlist(runlist)

    def get_data_store(self):
        """
        Return data store.

        """

        return self._data_store

    def get_observations(self, reflected_region=True):
        """
        Return observations.

        Reflected region analysis requires effective area and energy dispersion only.

        Parameters
        ----------
        reflected_region : bool
            Reflected region analysis.

        """

        available_irf = None
        if reflected_region:
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
            self._logger.error("Run list %s not found.", runlist)
            raise

        self._logger.info("Reading run list with %d observations from %s", len(_runs), runlist)
        return _runs

    def _from_target(self, target, viewcone):
        """
        Select data based on target coordinates and viewcone.

        Parameters
        ----------
        target : SkyCoord
            Target coordinates.
        viewcone : float
            Viewcone radius (deg).

        """

        observations = self._data_store.obs_table
        mask = target.separation(observations.pointing_radec) < viewcone * u.deg
        _runs = observations[mask]["OBS_ID"].data

        self._logger.info("Selecting %d runs from viewcone around %s", len(_runs), target)
        self._logger.info("WARNING - this is not tested")
        return _runs
