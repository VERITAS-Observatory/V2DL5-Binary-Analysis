"""
Orbital phase calculations

"""

import numpy as np


def get_orbital_phase(mjd, orbital_period=317.3, MJD0=54857.0, phase_reduce=True):
    """
    Calculate orbital phase for a given MJD.

    TODO: check how necessary the np.squeeze is

    Parameters
    ----------

    mjd: float
        MJD
    orbital_period: float
        Orbital period
    MJD0: float
        Reference MJD
    phase_reduce: bool
        Reduce phase to interval [0,1]

    Returns
    -------
    float
        Orbital phase

    """
    mjd = np.asarray(mjd)
    scalar_input = False
    if mjd.ndim == 0:
        mjd = mjd[None]  # Makes x 1D
        scalar_input = True

    phase = (mjd - MJD0) / orbital_period

    if phase_reduce:
        phase -= phase.astype(int)

    for i in range(len(phase)):
        if phase[i] < 0:
            phase[i] = 1.0 + phase[i]

    if scalar_input:
        return np.squeeze(phase)

    return phase


# def getMJDOrbitZeroPhase(mjd, object="HESS J0632+057", orbital_period=317.3):
#     """
#     Return MJD of phase 0 for the given orbital period period.
#
#     Parameters:
#     -----------
#     mjd: float
#         MJD
#     object: str
#         name of the object
#     orbital_period: float
#         orbital period
#
#     Returns:
#     --------
#     MJD0: float
#         MJD0 for the object
#
#     """
#     norbit = math.floor((mjd - get_mjd0()) / orbital_period)
#
#     return get_mjd0() + norbit * orbital_period
#
#
# def getNumberOfOrbits(mjd_min, mjd_max, Orbit=317.3):
#     """return number of orbits covered
#     Parameters:
#         - MJD (min, max value)
#         - orbital phase
#     Returns:
#         - number of orbits
#     """
#
#     return math.ceil((mjd_max - mjd_min) / Orbit)
#
#
#
# def getOrbitalPhaseRange(
#         mjd_min,
#         mjd_max,
#         phase_mean=0.5,
#         UppErr=True,
#         Orbit=317.3):
#     """calculate width of observation bin in orbital phase
#     Parameters:
#         - mjd_min, mjd_max: edges of MJD mean
#         - phase_mean: ??
#         - Up_err: upper error
#         - Orbit orbital phase
#     Returns:
#         - width of bin in orbital phase
#     """
#
#     ph_min = getOrbitalPhase(mjd_min, Orbit)
#     ph_max = getOrbitalPhase(mjd_max, Orbit)
#     if abs(ph_max - ph_min) < 1.e-3:
#         ph_err = 0.
#     elif ph_max > ph_min:
#         ph_err = 0.5 * (ph_max - ph_min)
#     else:
#         ph_err = 0.5 * (1. + ph_max - ph_min)
#
#     # check if error reaches over '0' or '1' -> simply cut (not entirely
#     # correct)
#     if UppErr:
#         if phase_mean + ph_err > 1.:
#             ph_err = 1. - phase_mean
#     else:
#         if phase_mean - ph_err < 0.:
#             ph_err = phase_mean
#
#     return ph_err
#
