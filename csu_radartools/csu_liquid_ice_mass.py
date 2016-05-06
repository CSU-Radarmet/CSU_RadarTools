"""
csu_liquid_ice_mass.py

Modifications by Timothy Lang
tjlangoc@gmail.com
1/27/2015
rev1 08/05/2015 - Python 3 compliant
rev2 11/20/2015 - Performance improvements

# Brody Fuchs, CSU, Oct 2014
# brfuchs@atmos.colostate.edu

# python remake of Brenda's code to calculate water and ice mass
# See original comments below
#***************************************************************************
#***************************************************************************
# This program will calculate ice and water mass at specified height levels
# given Z and ZDR at each grid point.
#
# Based on Code from A. Rowe and T. Lang from L. Nelson and L. Carey.
# Brenda Dolan
# February 27, 2012
# bdolan@atmos.colostate.edu
#***************************************************************************
#***************************************************************************

"""
from __future__ import division
import numpy as np
import warnings

DEFAULT_HFRZ = 4.5  # km MSL


def calc_liquid_ice_mass(dbz, zdr, z, T=None, Hfrz=None, method='cifelli',
                         fit_a=None, fit_b=None):
    """
    This function uses the ZDP method to discriminate between

    Arguments:

    Must be same shape
    ------------------
    dbz = Reflectivity (dBZ)
    zdr = Differential Reflectivity (dB)
    z = Height (km MSL)
    T = Temperature (deg C)

    Other keywords
    --------------
    method = Name of method to determine fit_a and fit_b for rain line.
             Options are 'cifelli' and 'cr1995' which correspond to
             Cifelli et al. 2002 and Carey and Rutledge (1995). If you put
             anything else here, it will use Carey and Rutledge (2000).
    fit_a, fit_b = Fit coefficients to the rain line. If user doesn't specify,
                   then the reference specified in method will be used.

    Scalar only
    -----------
    Hfrz = Height of freezing level (km MSL), if known; will be calculated
           from sounding provided otherwise.

    All need to be same array or scalar structure. If T == None, then
    function will assume default arrangement of delta_thresh.

    No error checking done, user is responsible for their own bad data masking.

    """
    len_flag = hasattr(dbz, '__len__')
    if len_flag:
        shape = np.shape(dbz)
        dbz = dbz.ravel()
        zdr = zdr.ravel()
        z = z.ravel()
        if T is not None:
            T = T.ravel()
    else:
        dbz = np.array([dbz])
        zdr = np.array([zdr])
        z = np.array([z])
        if T is not None:
            T = np.array([T])
    Hfrz = _get_hfrz(T, z, Hfrz)
    if Hfrz is None:
        return
    delta_thresh = get_delta_thresh(z, Hfrz)
    mass_w, mass_i = _liquid_ice_calcs(dbz, zdr, delta_thresh, fit_a=fit_a,
                                       fit_b=fit_b)
    if len_flag:
        return np.reshape(mass_w, shape), np.reshape(mass_i, shape)
    else:
        return mass_w[0], mass_i[0]


def linearize(dz):
    """dz = Reflectivity (dBZ), returns Z (mm^6 m^-3)"""
    return 10.0**(dz / 10.0)


def get_freezing_altitude(T, z):
    """
    T = temperature (deg C), z = altitude (any units)
    Assumes T & z are aligned, unmasked ndarrays of same shape.
    Returns highest possible freezing altitude.
    Partial soundings (i.e., no pass thru freezing altitude) get fitted to
    a regression line and Hfrz is then estimated via the intercept.
    """
    T = T.ravel()
    z = z.ravel()
    try:
        # Maybe there's already one or more zero points in sounding
        return np.max(z[T == 0])
    except:
        # Following gets indices of sorted z array, for comparison w/ T.
        # Commented two lines below were slower than necessary.
        # zs = sorted(z)
        # Tz = T[np.argsort(z)]
        sin = np.argsort(z)
        zs = z[sin]
        Tz = T[sin]

        # If sounding passes thru 0, then look immediately above and below
        if np.max(T) > 0 and np.min(T) < 0:
            zarg = np.argmax(zs[Tz > 0])
            T2 = Tz[zarg:zarg+2]
            z2 = zs[zarg:zarg+2]
            if T2[1] < T2[0]:
                return np.interp(0.0, T2[::-1], z2[::-1])
            else:
                return np.interp(0.0, T2, z2)
        # Sounding doesn't pass thru 0, trying linear regression
        elif np.max(T) < 0 or np.min(T) > 0:
            reg = np.polyfit(Tz, zs, 1)
            return reg[1]


def get_delta_thresh(z, Hfrz):
    """
    Makes assumptions about change in T with height,
    so an abnormal T profile may lead to inaccurate results.
    TJL - Currently feel set too aggressively, leading to turning obvious
    hail gates (high Z, low ZDR) into all-rain, giving huge LWCs.
    """
    delta_thresh = np.zeros_like(z) + 1.1
    cond = z <= Hfrz
    delta_thresh[cond] = (-8.0/3.0) * (z[cond] - Hfrz) + 2.0
    return delta_thresh


def get_linear_fits(method='cifelli'):
    """
    Get coefficients from linear fits to Zdp and Zh
    method = str, Reference for method
        'cifelli': Cifelli et al. (2002) - Amazon
        'cr1995': Carey and Rutledge (1995) - Colorado
        Default: Carey and Rutledge (2000) - Darwin
    """
    # TJL - Does Carey and Rutledge (1995) have another rain line?
    if method == 'cifelli':
        # These are Cifelli et al. (2002) for Amazon
        fit_a = 0.8178
        fit_b = 12.5088
    elif method == 'cr1995':
        # These are Carey and Rutledge (1995) for Colorado (?)
        fit_a = 0.9091
        fit_b = 8.5091
    else:
        # These are Carey and Rutledge (2000) for MCTEX
        fit_a = 0.77
        fit_b = 14.0
    return fit_a, fit_b


def calc_zh_zv(dbz, zdr):
    Zh = linearize(dbz)
    return Zh, Zh / linearize(zdr)  # Vertical component of reflectivity


def calc_zdp(Zh, Zv):
    Z_dp = Zh - Zv  # Difference reflectivity
    Z_dp[Z_dp <= 0] = 0.000000000001
    return 10.0 * np.log10(Z_dp)  # Convert to decibel units


def _get_hfrz(T, z, Hfrz):
    if T is None:
        if Hfrz is None:
            Hfrz = DEFAULT_HFRZ
    else:
        if Hfrz is None:
            if hasattr(T, '__len__'):
                Hfrz = get_freezing_altitude(T, z)
            else:
                warnings.warn('Given a scalar for T and no Hfrz, failing')
    return Hfrz


def _liquid_ice_calcs(dbz, zdr, delta_thresh, fit_a=None, fit_b=None):
    # Set up arrays to hold the output data
    mass_w = dbz * 0.0
    mass_i = dbz * 0.0

    # Get Zdp
    Zh, Zv = calc_zh_zv(dbz, zdr)
    Zdp = calc_zdp(Zh, Zv)

    # Determine ice fraction and partition reflectivity
    if fit_a is None:
        fit_a, fit_b = get_linear_fits(method='cifelli')
    dbz_est = fit_a * Zdp + fit_b  # Estimate the best fit line
    delta_z = dbz - dbz_est
    delta_z[delta_z < 0] = 0.0
    ice_frac = 1.0 - 10.0**(-0.1 * delta_z)
    Zh_rain = Zh * (1.0 - ice_frac)
    Zh_ice = Zh - Zh_rain

    # Check against delta_thresh
    check = delta_z - delta_thresh
    above = check >= 0
    below = check < 0

    # Z-M relationships for water and ice (Cifelli et al. 2002)
    mass_w[below] = 0.70 * 10.0**(-3) * Zh[below]**(0.886) * \
        (Zh[below] / Zv[below])**(-4.159)
    # mass_i[below] = 0.0 anyway
    mass_w[above] = 3.44 * 10.0**(-3) * Zh_rain[above]**(4.0/7.0)
    mass_i[above] = 1000.0 * np.pi * 917.0 * (4.0 * 10.0**6)**(3.0/7.0) * \
        ((5.28 * 10.0**(-18) * Zh_ice[above]) / 720.0)**(4.0/7.0)
    mass_w[zdr < 0] = 0.0  # Negative ZDR = all ice
    return mass_w, mass_i  # Units = g m^-3
