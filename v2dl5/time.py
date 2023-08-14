"""
Definition of time bins
"""

from astropy.time import Time


def get_list_of_nights(data_set, time_zone=-7.0):
    """
    Create a list of nights from a list of Time objects.

    Parameters
    ----------
    data_set : list
        List of Time objects.
    time_zone : float
        Time zone in hours (relative to UTC).
        Default: Mountain Standard Time (MST) = UTC-7

    Returns
    -------
    list
        List of nights (list of Time objects at midnight).

    """

    time_list = _get_starting_tims(data_set)

    # Reduce the list to intervals
    #    _mjd = set([(int)(time_obj.mjd) for time_obj in time_list])
    _mjd = {(int)(time_obj.mjd) for time_obj in time_list}

    time_intervals = []
    for _night in _mjd:
        _mid_night = float(_night) - time_zone / 24.0
        time_intervals.append(Time([_mid_night - 0.5, _mid_night + 0.5], format="mjd", scale="utc"))

    return time_intervals


def _get_starting_tims(data_set):
    """
    Get the starting times of the observations.

    Parameters
    ----------
    data_set : `~gammapy.data.DataStoreObservation`
        Data set.

    Returns
    -------
    list
        List of starting times.

    """

    _time_start = []
    for _data_set in data_set:
        _time_start.extend(_data_set.gti.time_start)

    return _time_start
