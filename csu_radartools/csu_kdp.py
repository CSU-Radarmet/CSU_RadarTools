"""
Timothy James Lang
tjlangco@gmail.com

Last Updated 03 May 2016 (Python 2.7/3.4)
Last Updated 26 July 2005 (IDL)

csu_kdp v1.6

Change Log
----------
v1.6 Major Changes (05/03/2016):
1. Cython now an option for speeding up KDP routines.

v1.5.1 Major Changes (12/08/2015):
1. Made number of gates used for standard deviation calculation adjustable.

v1.5 Major Changes (11/17/2015):
1. Now using a Fortran shared object (calc_kdp_ray_fir) to
   do the ray-based KDP calculations. This has vastly sped up
   the overall KDP processing (> 100x). f2py FTW!
2. Updated calc_kdp_bringi to fail softly when (window/gs) is not
   even.

v1.4 Major Changes (09/04/2015):
1. Added window keyword to enable stretching the FIR window (e.g.,
   use a 21-pt filter over 5 km with 250-m gate spacing).
2. Forcing FIR order to be even, _calc_kdp_ray will crash otherwise.

v1.3 Major Changes (08/05/2015):
1. Made Python 3 compatible.
2. Fixed issue with non-integer array indices.

v1.2 Major Changes (07/10/2015):
1. Made sub-module pep8 compliant.

v1.1 Major Changes (04/27/2015):
1. Made algorithm work with a user-defined gate spacing (via gs keyword).
   Untested on gate spacings that do not divide evenly into the 3-km window
   used for filtering the PHIDP data, however. But common gate spacings
   like 50, 100, 150, 200, 250, and 300 meters should all work fine.
2. Made the algorithm capable of receiving 2D array inputs (i.e., azimuth &
   range) as well as 1D inputs (range only). If 2D, rng needs to be 2D as
   well. However, thsd should remain a scalar, or 1D and only vary by range.

"""
from __future__ import division, print_function
import numpy as np
from numpy import linalg
from scipy.signal import firwin
from warnings import warn
from calc_kdp_ray_fir import calc_kdp_ray_fir
# import time

VERSION = '1.6'

# Used by FIR coefficient function (get_fir)
FIR_GS = 150.0
FIR_WIN = 3.0
FIR_ORDER = 20
FIR_GAIN = 1.0
FIR_FREQ = 0.08
FIR_STD = 28.0
KM2M = 1000.0
STD_GATE = 11


def calc_kdp_bringi(dp=None, dz=None, rng=None, thsd=12, nfilter=1,
                    bad=-32768, gs=FIR_GS, window=FIR_WIN, std_gate=STD_GATE):
    """
    Overview
    --------
    This is an old algorithm that uses an FIR filter to process differential
    phase and extract specific differential phase. It works on polarimetric
    radar data. It is based on code provided by V. N. Bringi and Yanting Wang
    of CSU Electrical Engineering. It assumes differential phase has been
    unfolded already. You can send this function either 1D or 2D arrays of
    data. If 2D, it assumes the first index is azimuth so it will loop over
    that, calculating KDP along individual rays.

    Steps
    -----
    1. Standard deviation of differential phase is calculated and used to
       QC the phase data. The stdev calculation uses up to std_gate consecutive
       gates regardless of gate spacing.
    2. Differential phase is filtered using the FIR filter, which has been
       tuned to the number of gates contained within the FIR window. This
       algorithm only works for window / gate spacing = even number.
    3. Specific differential phase is calculated by consulting reflectivity.
       As reflectivity declines progressively more and more gates are needed
       in the window used to fit a line to the filtered phase. Specific
       differential phase is half the slope of that line.

    Reference
    ---------
    Timothy J. Lang, David A. Ahijevych, Stephen W. Nesbitt, Richard E.
    Carbone, Steven A. Rutledge, and Robert Cifelli, 2007: Radar-Observed
    Characteristics of Precipitating Systems during NAME 2004. J. Climate,
    20, 1713–1733. doi: http://dx.doi.org/10.1175/JCLI4082.1

    Arguments
    ---------
    dp = Differential phase (deg, 1D or 2D array)
    dz = Reflectivity (dBZ, 1D or 2D array)
    rng = Range (km, 1D or 2D array -
          use np.meshgrid() first tp make rng 2D if needed)
    thsd = Threshold for standard deviation of differential phase, above which
           the data are not considered when filtering or calculating specific
           differential phase. The user can specify a 1D vector of spatially
           varying thresholds instead (i.e., vary by range).
    nfilter = Number of times to apply the FIR filter
    bad = Value for bad/missing data
    gs = Gate spacing of radar (meters)
    window = Changes window over which FIR filter is applied (km). Also affects
             the width of the adaptive KDP calculations.
    std_gate = Number of gates for standard deviation of phase calculation.
               Must be odd or function will just set it to the default value.

    Returns
    -------
    kd_lin = Specific differential phase (deg/km, 1D or 2D array)
    dp_lin = Filtered differential phase (deg, 1D or 2D array)
    sd_lin = Standard deviation of diff. phase (deg, 1D or 2D array)

    """
    # Quick check on all vars. Used keywords so order doesn't matter.
    if dp is None or dz is None or rng is None:
        warn('Missing needed variables (dp, dz, and/or rng), failing ...')
        return
    if np.ndim(dp) != np.ndim(dz) or np.ndim(dp) != np.ndim(rng):
        warn('Array sizes don\'t match, failing ...')
        return
    if std_gate % 2 != 1:
        warn('std_gate must be odd, using ' + str(STD_GATE) +
             ' gates as the default window')
        std_gate = STD_GATE
    fir = get_fir(gs=gs, window=window)
    if fir is None:
        print('Fix window/gs to be even, failing ...')
        return None, None, None

    # Following lines ensure right dtype passed to Cython extensions (if used)
    dp = np.array(dp).astype('float32')
    dz = np.array(dz).astype('float32')
    rng = np.array(rng).astype('float32')
    fir['coef'] = np.array(fir['coef']).astype('float64')

    if not hasattr(thsd, '__len__'):
        thsd = np.zeros_like(dp) + thsd
    # If array is 2D, then it assumes the first index refers to azimuth.
    # Thus it loops over that.
    if np.ndim(dp) == 2:
        kd_lin = np.zeros_like(dp) + bad
        dp_lin = np.zeros_like(dp) + bad
        sd_lin = np.zeros_like(dp) + 100.0
        for ray in np.arange(np.shape(dp)[0]):
            dpl = len(dp[ray])
            kd_lin[ray], dp_lin[ray], sd_lin[ray] = calc_kdp_ray_fir(
                dpl, dp[ray], dz[ray], rng[ray], thsd[ray],
                nfilter, bad, fir['order'], fir['gain'], fir['coef'], std_gate)
#             kd_lin[ray], dp_lin[ray], sd_lin[ray] = \
#                 _calc_kdp_ray(dp[ray], dz[ray], rng[ray], thsd=thsd,
#                               nfilter=nfilter, bad=bad, fir=fir)
    # Or
    elif np.ndim(dp) == 1:
        kd_lin, dp_lin, sd_lin = calc_kdp_ray_fir(
            len(dp), dp, dz, rng, thsd, nfilter, bad,
            fir['order'], fir['gain'], fir['coef'], std_gate)
#        kd_lin, dp_lin, sd_lin = _calc_kdp_ray(
#             dp, dz, rng, thsd=thsd, fir=fir, nfilter=nfilter, bad=bad)
    else:
        warn('Need 2D or 1D array, failing ...')
        return
    return kd_lin, dp_lin, sd_lin


def get_fir(gs=FIR_GS, window=FIR_WIN):
    """
    gs = Gate Spacing (m)
    window = Filter Window (km)
    window divided by gs should be an even number!
    """
    fir = {}
    fir['order'] = np.int32(window * KM2M / gs)
    if fir['order'] % 2 != 0:
        warn('gs / window must be an even number! #Failing ...')
        return
    fir['gain'] = FIR_GAIN
    # ratio = FIR_GS / gs
    ratio = fir['order'] / FIR_ORDER
    freq = FIR_FREQ / ratio
    std = ratio * FIR_STD
    fir['coef'] = firwin(fir['order'] + 1, freq, window=('gaussian', std))
    # print('debug', fir)
    return fir


def _calc_kdp_ray(dp, dz, rng, thsd=12, nfilter=1, bad=-32768, fir=None):
    """
    Pure Python approach to filtering phase and estimating KDP. Currently
    disabled due to performance issues.

    Arguments
    ---------
    dp = 1D ray of differential phase
    dz = 1D ray of reflectivity
    rng = 1D ray of range
    thsd = Scalar or 1D ray of diff phase standard deviation thresholds
    nfilter = Number of times to filter the data
    bad = Bad/missing data value
    fir = Dictionary containing FIR filter parameters

    Returns
    -------
    kd_lin = Specific differential phase (deg/km, 1D array)
    dp_lin = Filtered differential phase (deg, 1D array)
    sd_lin = Standard deviation of diff. phase (deg, 1D array)
    """
    # Define needed variables
    kd_lin = np.zeros_like(rng) + bad
    sd_lin = np.zeros_like(rng) + 100.0
    # User can provide a spatially varying stddev(dp) threshold
    if not hasattr(thsd, '__len__'):
        thsd = np.zeros_like(rng) + thsd
    length = len(rng)
    lin = np.arange(length)
    # Half window size for calculating stdev of phase (fixed @ 11 gates)
    half_std_win = 5
    half_fir_win = fir['order'] // 2  # Half window size for FIR filtering
    y = np.zeros(length) + bad  # Dummy variable to store filtered phase
    z = 1.0 * dp  # Dummy variable to store un/pre-processed phase
    # print(time.time() - begin_time, 'seconds since start (DEF)')

    #####################################################################
    # Calculate standard deviation of phidp
    mask = dp != bad
    for i in lin[mask]:
        index1 = np.int32(i - half_std_win)
        index2 = np.int32(i + half_std_win)
        if index1 >= 0 and index2 < length - 1:
            yy = dp[index1:index2]
            tmp_mask = mask[index1:index2]
            if len(yy[tmp_mask]) > half_std_win:
                sd_lin[i] = _quick_std(yy, tmp_mask)

    # ------------- MAIN LOOP of Phidp Adaptive Filtering ------------------
    # FIR FILTER SECTION
    for mloop in np.arange(nfilter):
        mask = np.logical_and(sd_lin <= thsd, z != bad)
        for i in lin[mask]:
            index1 = np.int32(i - half_fir_win)
            index2 = np.int32(i + half_fir_win)
            if index1 >= 0 and index2 < length - 1:
                yy = z[index1:index2+1]
                xx = rng[index1:index2+1]
                tmp_mask = mask[index1:index2+1]
                siz = len(yy[tmp_mask])
                if siz > 0.8 * fir['order']:
                    if siz < fir['order'] + 1:
                        result = _leastsqrs(xx, yy, siz, tmp_mask)
                        yy[~tmp_mask] = result[0] * xx[~tmp_mask] + result[1]
                    y[i] = fir['gain'] * np.dot(fir['coef'], yy)
        z = 1.0 * y  # Enables re-filtering of processed phase
    dp_lin = 1.0 * z
    # print(time.time() - begin_time, 'seconds since start (FDP)')
    # *****************END LOOP for Phidp Adaptive Filtering******************

    # CALCULATE KDP
    # Default value for nadp is half_fir_win, but varies based on Zh
    nadp = np.int16(0 * dz + half_fir_win)
    tmp_mask = dz < 35
    nadp[tmp_mask] = 3 * half_fir_win
    tmp_mask = np.logical_and(dz >= 35, dz < 45)
    nadp[tmp_mask] = 2 * half_fir_win
    mask = dp_lin != bad
    for i in lin[mask]:
        index1, index2 = _get_nadp_indices(nadp, i)
        if index1 >= 0 and index2 <= length:
            tmp_mask = mask[index1:index2]
            xx = rng[index1:index2]
            siz = len(xx[tmp_mask])
            # Improved Kdp based on LSE fit to Adap filt Phidp
            if siz >= 0.8 * nadp[i]:
                yy = dp_lin[index1:index2]
                kd_lin[i] = _fit_line_and_get_kdp(xx, yy, siz, tmp_mask)
    # *******************END KDP CALCULATION****************************
    # print(time.time() - begin_time, 'seconds since start (KDP/Done)')
    return kd_lin, dp_lin, sd_lin


def _leastsqrs(xx, yy, siz, tmp_mask):
    """
    Following is faster than np.polyfit
    e.g., return np.polyfit(xx[tmp_mask], yy[tmp_mask], 1)
    """
    A = np.array([xx[tmp_mask], np.ones(siz)])
    return linalg.lstsq(A.T, yy[tmp_mask])[0]


def _get_nadp_indices(nadp, i):
    half_nadp = nadp[i] / 2
    return np.int32(i - half_nadp), np.int32(i + half_nadp + 1)


def _fit_line_and_get_kdp(xx, yy, siz, tmp_mask):
    result = _leastsqrs(xx, yy, siz, tmp_mask)
    return 0.5 * result[0]


def _quick_std(array, mask):
    """Following is faster than np.std()"""
    a = array[mask]
    m = a.mean()
    c = a - m
    return (np.dot(c, c) / a.size)**0.5
