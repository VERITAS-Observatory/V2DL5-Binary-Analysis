"""
Data class holding data store and observations.
"""

import logging

import numpy as np
from astropy import units as u
from gammapy.data import DataStore

import v2dl5.bti as BTI


class Data:
    """
    Data class holding data store and observations

    Allows to select data from run list or based on
    target coordinates (and observation cone).

    Parameters
    ----------
    run_list : str
        Path to run list.
    data_directory : str
        Path to data directory (holding hdu-index.fits.gz and obs-index.fits.gz).
    target : SkyCoord
        Target coordinates.
    obs_cone_radius : float
        observation cone radius (deg).

    """

    def __init__(self, args_dict, target=None):
        """
        Initialize Data object.

        Uses run_list if not set to None, otherwise selects data
        according to target coordinates and observation cone.

        """

        self._logger = logging.getLogger(__name__)

        self._logger.info(
            "Initializing data object from %s", args_dict["observations"]["datastore"]
        )
        self._data_store = DataStore.from_dir(args_dict["observations"]["datastore"])
        self.target = target
        if args_dict.get("run_list") is None:
            self.runs = self._from_target(
                args_dict["observations"].get("obs_cone_radius", 5.0 * u.deg)
            )
        else:
            self.runs = self._from_run_list(args_dict.get("run_list"))
        self._update_gti(args_dict.get("bti", None))

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

        required_irf = "point-like"
        return self._data_store.get_observations(self.runs, required_irf=required_irf)

    def _from_run_list(self, run_list):
        """
        Read run_list from file and select data

        Parameters
        ----------
        run_list : str
            Path to run list.

        """

        if run_list is None:
            return None

        try:
            _runs = np.loadtxt(run_list, dtype=int)
        except OSError:
            self._logger.error("Run list %s not found.", run_list)
            raise

        self._logger.info("Reading run list with %d observations from %s", len(_runs), run_list)
        return _runs

    def _from_target(self, obs_cone_radius):
        """
        Select data based on target coordinates and observation cone.

        Parameters
        ----------
        obs_cone_radius : float
            observation cone radius (deg).

        """

        observations = self._data_store.obs_table
        mask = self.target.separation(observations.pointing_radec) < obs_cone_radius * u.deg
        _runs = observations[mask]["OBS_ID"].data

        self._logger.info(
            "Selecting %d runs from observation cone around %s", len(_runs), self.target
        )
        self._logger.warning("THIS IS NOT TESTED")
        return _runs

    def get_on_region_radius(self):
        """
        Return on region radius.

        Simplest case. Ignores possible energy and offset dependence.

        """

        observations = self.get_observations()
        try:
            rad_max = set(obs.rad_max.data[0][0] for obs in observations)
        except IndexError:
            self._logger.error("Rad max not found in observations.")
            raise

        if len(rad_max) > 1:
            self._logger.error("Rad max is not the same for all observations.")
            raise ValueError

        on_region = rad_max.pop() * u.deg
        self._logger.info(f"On region size: {on_region}")

        return on_region

    def get_max_wobble_distance(self, fov=3.5 * u.deg):
        """
        Return maximum distance from target position.
        Add if necessary the telescope field of view.

        Parameters
        ----------
        fov : astropy.units.Quantity
            Telescope field of view.

        Returns
        -------
        max_offset : astropy.units.Quantity
            Maximum offset.

        """

        woff = np.array(
            [
                self.target.separation(obs.pointing.get_icrs()).degree
                for obs in self.get_observations()
            ]
        )

        return np.max(woff) * u.deg + fov / 2.0

    def _update_gti(self, bti):
        """
        Update good time intervals by removing bad time intervals.

        Parameters
        ----------
        bti : list of dict
            List of bad time intervals
            Given us {"run": run, "bti_start": start, "bti_length": length}

        """

        if bti is None:
            return

        for obs in self.get_observations():
            bti_pairs = [
                (item["bti_start"], item["bti_start"] + item["bti_length"])
                for item in bti
                if item["run"] == obs.obs_id
            ]
            if len(bti_pairs) == 0:
                self._logger.debug("No BTI found for {obs.obs_id}")
                continue
            self._logger.debug("Updating GTI for {obs.obs_id} with {bti_pairs}")
            obs.gti.stack(other=BTI.BTI(obs).update_gti(bti_pairs))
            obs.gti.union()
