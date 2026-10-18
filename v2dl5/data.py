"""Data class holding data store and observations."""

import logging
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.table import Table
from gammapy.data import DataStore

import v2dl5.bti as BTI  # noqa: N812


class Data:
    """
    Data class holding data store and observations.

    Allows to select data from run list or based on
    target coordinates and observation cone.

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

        Uses 'run_list' from args_dict if not set to None, otherwise selects data
        according to target coordinates and observation cone.

        """
        self._logger = logging.getLogger(__name__)
        self._observation_cache = {}
        self._bti = args_dict.get("bti")

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

    def preflight(self, output_file):
        """Read required HDUs and account for every requested run before reduction."""
        rows = []
        observations = []
        for run in self.runs:
            try:
                if np.count_nonzero(self._data_store.obs_table["OBS_ID"] == run) != 1:
                    raise ValueError("OBS_ID is missing or ambiguous in observation index")
                obs = self._data_store.get_observations(
                    [run], required_irf="point-like", skip_missing=False
                )[0]
                selected = [obs]
                self._apply_bti(selected)
                obs = selected[0]
                events = obs.events
                if obs.gti.time_sum <= 0 * u.s:
                    raise ValueError("No usable GTI")
                if np.any(obs.gti.time_stop <= obs.gti.time_start) or np.any(
                    obs.gti.time_start[1:] < obs.gti.time_stop[:-1]
                ):
                    raise ValueError("Invalid or overlapping GTIs")
                in_gti = np.zeros(len(events.table), dtype=bool)
                for low, high in obs.gti.time_intervals:
                    in_gti |= (events.time >= low) & (events.time < high)
                if not np.all(in_gti):
                    raise ValueError("Events outside GTIs; reconcile preprocessing and time cuts")
                if np.any(events.energy <= 0 * u.TeV):
                    raise ValueError("Nonpositive event energy")
                livetime = obs.observation_live_time_duration.to_value("s")
                if not np.isfinite(livetime) or livetime <= 0:
                    raise ValueError("Invalid livetime")
                if not obs.aeff.is_pointlike or not obs.edisp.is_pointlike:
                    raise ValueError("Point-like IRFs required")
                for irf in (obs.aeff, obs.edisp):
                    if not np.all(np.isfinite(irf.data)) or np.any(irf.data < 0):
                        raise ValueError("Invalid IRF values")
                if not np.any(obs.aeff.data > 0) or not np.any(obs.edisp.data > 0):
                    raise ValueError("Empty IRF")
                if obs.rad_max is None:
                    raise ValueError("Missing RAD_MAX")
                offset = self.target.separation(obs.pointing.get_icrs()).to_value("deg")
                bounds = obs.aeff.axes["offset"].edges.to_value("deg")
                if not np.isfinite(offset) or not bounds[0] <= offset <= bounds[-1]:
                    raise ValueError("Target offset outside effective-area support")
                for name in ("RA", "DEC", "ENERGY"):
                    if not np.all(np.isfinite(events.table[name])):
                        raise ValueError(f"Nonfinite events: {name}")
                if not np.all(np.isfinite(events.time.jd)):
                    raise ValueError("Nonfinite event times")
                observations.append(obs)
                rows.append((int(run), "available", "", livetime, offset))
            except (OSError, ValueError, KeyError, IndexError, TypeError) as exc:
                rows.append((int(run), "missing_data", str(exc), 0.0, float("nan")))
        Table(rows=rows, names=("obs_id", "status", "reason", "livetime_s", "offset_deg")).write(
            output_file, format="ascii.ecsv", overwrite=True
        )
        if len(observations) != len(self.runs):
            raise ValueError("DL3 preflight failed; see availability.ecsv")
        self._observation_cache[(True, False)] = observations

    def get_data_store(self):
        """Return data store."""
        return self._data_store

    def get_observations(self, reflected_region=True, skip_missing=False):
        """
        Return list of observations.

        Parameters
        ----------
        reflected_region : bool
            Reflected region analysis.
        skip_missing : bool
            Skip missing observations.

        Returns
        -------
        observations : list of `~gammapy.data.Observation`
            List of observations.

        """
        cache_key = (reflected_region, skip_missing)
        if cache_key not in self._observation_cache:
            required_irf = "point-like" if reflected_region else "full-enclosure"
            self._observation_cache[cache_key] = self._data_store.get_observations(
                self.runs,
                required_irf=required_irf,
                skip_missing=skip_missing,
            )
            self._apply_bti(self._observation_cache[cache_key])
        return self._observation_cache[cache_key]

    @staticmethod
    def _from_run_list(run_list):
        """
        Read run list from file and select data.

        Parameters
        ----------
        run_list : str
            Path to run list.

        """
        if run_list is None:
            return None

        runs = []
        for number, line in enumerate(Path(run_list).read_text().splitlines(), 1):
            value = line.split("#", 1)[0].strip()
            if not value:
                continue
            if not value.isdecimal() or int(value) <= 0:
                raise ValueError(f"Invalid run ID at {run_list}:{number}: {value}")
            run = int(value)
            if run in runs:
                raise ValueError(f"Duplicate run {run} in {run_list}")
            runs.append(run)
        if not runs:
            raise ValueError(f"Run list is empty: {run_list}")
        return runs

    def _from_target(self, obs_cone_radius):
        """
        Select data based on target coordinates and observation cone.

        Parameters
        ----------
        obs_cone_radius : float
            observation cone radius (deg).

        """
        if self.target is None:
            raise ValueError("A target is required for target-based data selection")
        observations = self._data_store.obs_table
        obs_cone_radius = u.Quantity(obs_cone_radius).to(u.deg)
        mask = self.target.separation(observations.pointing_radec) < obs_cone_radius
        _runs = observations[mask]["OBS_ID"].data

        self._logger.info(
            "Selecting %d runs from observation cone around %s", len(_runs), self.target
        )
        if len(_runs) == 0:
            raise ValueError(f"No observations found within {obs_cone_radius} of {self.target}")
        return _runs

    def get_on_region_radius(self):
        """
        Return on region radius.

        Simplest case. Ignores possible energy and offset dependence.

        """
        radii = []
        for obs in self.get_observations():
            if obs.rad_max is None:
                raise ValueError("Missing RAD_MAX")
            values = obs.rad_max.quantity.to_value(u.deg)
            if not np.all(np.isfinite(values)) or np.any(values <= 0):
                raise ValueError("Invalid RAD_MAX")
            if not np.allclose(values, values.flat[0], rtol=1e-6, atol=0):
                raise ValueError("Only constant RAD_MAX is supported")
            radii.append(values.flat[0])
        if not radii or not np.allclose(radii, radii[0], rtol=1e-6, atol=0):
            raise ValueError("A common constant RAD_MAX is required; split radius groups")
        return radii[0] * u.deg

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
            Maximum offset (radius of FoV).

        """
        woff = np.array(
            [
                self.target.separation(obs.pointing.get_icrs()).degree
                for obs in self.get_observations()
            ]
        )
        return np.max(woff) * u.deg + u.Quantity(fov) / 2.0

    def _update_gti(self, bti):
        """
        Update good time intervals by removing bad time intervals.

        Parameters
        ----------
        bti : list of dict
            List of bad time intervals
            Given us {"run": run, "bti_start": start, "bti_length": length}

        """
        if bti is not None:
            self._bti = bti

        for observations in self._observation_cache.values():
            self._apply_bti(observations)

    def _apply_bti(self, observations):
        """Apply configured bad-time intervals to an observation collection."""
        if self._bti is None:
            return

        for index, obs in enumerate(observations):
            bti_pairs = [
                (item["bti_start"], item["bti_start"] + item["bti_length"])
                for item in self._bti
                if item["run"] == obs.obs_id
            ]
            if len(bti_pairs) == 0:
                self._logger.debug("No BTI found for %s", obs.obs_id)
                continue
            self._logger.debug("Updating GTI for %s with %s", obs.obs_id, bti_pairs)
            gti = BTI.BTI(obs).update_gti(bti_pairs)
            events = getattr(obs, "events", None)
            if events is None:
                obs.gti = gti
            else:
                keep = np.zeros(len(events.table), dtype=bool)
                for low, high in gti.time_intervals:
                    keep |= (events.time >= low) & (events.time < high)
                observations[index] = obs.copy(
                    in_memory=True, gti=gti, events=events.select_row_subset(keep)
                )
