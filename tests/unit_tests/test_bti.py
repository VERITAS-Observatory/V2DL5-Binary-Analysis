import logging

import astropy.units as u
import pytest
from astropy.time import Time
from gammapy.data import GTI, Observation

import v2dl5.bti as BTI

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


@pytest.fixture
def get_good_gti():
    gti = GTI.create(
        [0 * u.s],
        [80 * u.s],
        reference_time=Time("2023-10-15T16:04"),
    )
    return gti


@pytest.fixture
def get_gti():
    gti = GTI.create(
        [0 * u.s, 10 * u.s, 20 * u.s],
        [5 * u.s, 15 * u.s, 25 * u.s],
        reference_time=Time("2023-10-15T16:04"),
    )
    return gti


class TestMaxTimeInterval:
    # Return the duration between the beginning of the first and the end of the last GTI.
    def test_duration_between_first_and_last_GTI(self, get_gti):
        obs = Observation(gti=get_gti)
        bti = BTI.BTI(obs=obs)
        duration = bti.max_time_interval()
        assert duration == 25

    # Return an integer value.
    def test_return_integer_value(self, get_gti):
        obs = Observation(gti=get_gti)
        bti = BTI.BTI(obs=obs)
        duration = bti.max_time_interval()
        assert isinstance(duration, int)

    # The GTI object has only one GTI.
    def test_one_GTI(self):
        obs = Observation(gti=GTI.create([0 * u.s], [5 * u.s], reference_time=None))
        bti = BTI.BTI(obs=obs)
        duration = bti.max_time_interval()
        assert duration == 5

    # The GTI object has GTIs with the same start and stop times.
    def test_same_start_and_stop_times(self):
        obs = Observation(
            gti=GTI.create(
                [0 * u.s, 10 * u.s, 20 * u.s], [0 * u.s, 10 * u.s, 20 * u.s], reference_time=None
            )
        )
        bti = BTI.BTI(obs=obs)
        duration = bti.max_time_interval()
        assert duration == 20


class Test_CreateOnTimeList:
    # Create a list of on and off times in 1 s bin with no bad time intervals.
    def test_no_bad_time_intervals(self, get_good_gti):
        obs = Observation(gti=get_good_gti)
        bti = []
        bti_obj = BTI.BTI(obs)
        on_time_list = bti_obj._create_on_time_list(bti)
        assert len(on_time_list) == bti_obj.max_time_interval() + 1
        assert all(x == 1 for x in on_time_list)

    # Create a list of on and off times in 1 s bin with one bad time interval.
    def test_one_bad_time_interval(self, get_good_gti):
        obs = Observation(gti=get_good_gti)
        bti = [(10, 20)]
        bti_obj = BTI.BTI(obs)
        on_time_list = bti_obj._create_on_time_list(bti)
        assert len(on_time_list) == bti_obj.max_time_interval() + 1
        assert on_time_list[10:21] == [0] * 11
        assert all(x == 1 for x in on_time_list[:10] + on_time_list[21:])

    # Create a list of on and off times in 1 s bin with multiple non-overlapping bad time intervals.
    def test_multiple_non_overlapping_bad_time_intervals(self, get_good_gti):
        obs = Observation(gti=get_good_gti)
        bti = [(10, 20), (30, 40), (50, 60)]
        bti_obj = BTI.BTI(obs=obs)
        on_time_list = bti_obj._create_on_time_list(bti)
        assert len(on_time_list) == bti_obj.max_time_interval() + 1
        assert on_time_list[10:21] == [0] * 11
        assert on_time_list[30:41] == [0] * 11
        assert on_time_list[50:61] == [0] * 11
        assert all(
            x == 1
            for x in on_time_list[:10]
            + on_time_list[21:30]
            + on_time_list[41:50]
            + on_time_list[61:]
        )

    # Create a list of on and off times in 1 s bin with multiple overlapping bad time intervals.
    def test_multiple_overlapping_bad_time_intervals(self, get_good_gti):
        obs = Observation(gti=get_good_gti)
        bti = [(10, 20), (15, 40), (50, 60)]
        bti_obj = BTI.BTI(obs=obs)
        on_time_list = bti_obj._create_on_time_list(bti)
        assert len(on_time_list) == bti_obj.max_time_interval() + 1
        assert on_time_list[10:41] == [0] * 31
        assert all(x == 1 for x in on_time_list[:10] + on_time_list[41:50] + on_time_list[61:])

    # Create a list of on and off times in 1 s bin with bad time intervals that start before
    # the observation start time.
    def test_bad_time_intervals_start_before_observation_start(self, get_good_gti):
        obs = Observation(gti=get_good_gti)
        bti = [(-10, 10), (70, 90), (100, 110), (-100, -50)]
        bti_obj = BTI.BTI(obs)
        on_time_list = bti_obj._create_on_time_list(bti)
        assert len(on_time_list) == bti_obj.max_time_interval() + 1
        assert on_time_list[:11] == [0] * 11
        assert all(x == 1 for x in on_time_list[11:69])


class Test_ExtractGTIStartStopTimes:
    # Create a list of on and off times in 1 s bin with no bad time intervals.
    def test_no_bad_time_intervals(self, get_good_gti):
        obs = Observation(gti=get_good_gti)
        bti = []
        bti_obj = BTI.BTI(obs)
        on_time_list = bti_obj._create_on_time_list(bti)
        start, stop = bti_obj._extract_gti_start_stop(on_time_list)
        assert start == 0.0 * u.s
        assert stop == 80 * u.s

    # Create a list of on and off times in 1 s bin with multiple overlapping bad time intervals.
    def test_multiple_overlapping_bad_time_intervals(self, get_good_gti):
        obs = Observation(gti=get_good_gti)
        bti = [(10, 20), (15, 40), (50, 60)]
        bti_obj = BTI.BTI(obs=obs)
        on_time_list = bti_obj._create_on_time_list(bti)
        start, stop = bti_obj._extract_gti_start_stop(on_time_list)
        assert start == [0.0 * u.s, 41.0 * u.s, 61.0 * u.s]
        assert stop == [9 * u.s, 49 * u.s, 80 * u.s]


class Test_MetTstart:
    # Return start time in seconds since MET.
    def test_return_start_time(self, get_good_gti):
        obs = Observation(gti=get_good_gti)
        bti = BTI.BTI(obs)
        result = bti._met_tstart()
        assert isinstance(result, u.Quantity)
        assert result.value == 0.0
