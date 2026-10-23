"""Regression checks for selecting events within good-time intervals."""

import copy
import logging
from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from gammapy.data import GTI, EventList

from v2dl5.data import Data


@pytest.mark.parametrize("times", [
    [-1, 0, 1, 2, 3, 4, 5, 5.023883],
    [-1, 2, 3, 5, 5.023883],
])
def test_preflight_filters_events_without_rejecting_run(tmp_path, times):
    reference = Time("2023-01-01", scale="tt")
    gti = GTI.create([0, 3] * u.s, [2, 5] * u.s, reference_time=reference)
    events = EventList(Table({
        "TIME": reference + np.array(times) * u.s,
        "RA": np.zeros(len(times)) * u.deg,
        "DEC": np.zeros(len(times)) * u.deg,
        "ENERGY": np.ones(len(times)) * u.TeV,
    }))

    class Observation(SimpleNamespace):
        def copy(self, in_memory=False, **kwargs):
            result = copy.copy(self)
            result.__dict__.update(kwargs)
            return result

    irf = SimpleNamespace(is_pointlike=True, data=np.ones(1), axes={
        "offset": SimpleNamespace(edges=[0, 1] * u.deg),
    })
    obs = Observation(
        events=events, gti=gti, aeff=irf, edisp=irf, rad_max=object(),
        observation_live_time_duration=3.6 * u.s,
        pointing=SimpleNamespace(get_icrs=lambda: SkyCoord(0, 0, unit="deg")),
    )
    data = Data.__new__(Data)
    data.runs = [7]
    data.target = SkyCoord(0, 0, unit="deg")
    data._logger = logging.getLogger(__name__)
    data._bti = None
    data._observation_cache = {}
    data._data_store = SimpleNamespace(
        obs_table=Table({"OBS_ID": [7]}),
        get_observations=lambda *_args, **_kwargs: [obs],
    )
    output = tmp_path / "availability.ecsv"
    data.preflight(output)
    selected = data.get_observations()[0]
    expected = [t for t in times if 0 <= t < 2 or 3 <= t < 5]
    np.testing.assert_allclose(
        (selected.events.time - reference).to_value("s"), expected, atol=1e-6,
    )
    assert len(obs.events.table) == len(times)
    assert selected.gti.time_sum.to_value("s") == pytest.approx(4)
    assert selected.observation_live_time_duration.to_value("s") == pytest.approx(3.6)
    row = Table.read(output)[0]
    assert row["status"] == "available"
    assert row["livetime_s"] == pytest.approx(3.6)
