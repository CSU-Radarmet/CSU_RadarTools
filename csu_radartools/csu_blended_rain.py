"""
csu_blended_rain.py

# Brody Fuchs, CSU, Oct 2014
# brfuchs@atmos.colostate.edu

# python remake of Brenda's code to calculate water and ice mass

Amendments by
Timothy Lang (tjlangco@gmail.com)
2/20/2015
rev1 07/10/2015
rev2 08/03/2015 - Python 3
rev3 09/16/2015 - Fixed logical inconsistencies leading to lack of
                  rainfall calculation in HID = rain + low Z + high Kdp/Zdr
rev4 05/10/2016 - Added ability for user to provide custom parameters to
                  polarimetric rainfall equations via the blended rainfall
                  routines. Also customized blended routines to handle
                  non-S bands. In this case, only R-Z and R-Kdp are used.
"""

from __future__ import absolute_import
from __future__ import division
import numpy as np
import warnings
from .csu_liquid_ice_mass import (calc_zh_zv, calc_zdp,
                                  linearize, get_linear_fits)


def calc_rain_zr(dz, a=300.0, b=1.4):
    """
    dz = Reflectivity (dBZ), returns rainfall rate (mm h^-1)
    Form Z = a * R**b
    """
    return (linearize(dz) / a)**(1.0 / b)


def calc_rain_nexrad(dz, thresh_nexrad=53.0, a=300.0, b=1.4):
    """
    dz = Reflectivity (dBZ)
    thresh_nexrad = Reflectivity cap for rain calc
    a, b => Z = a * R**b
    """
    rr = calc_rain_zr(dz, a=a, b=b)
    cond = dz > thresh_nexrad
    rr[cond] = calc_rain_zr(thresh_nexrad, a=a, b=b)
    return rr


def calc_rain_kdp_zdr(kdp, zdr, a=90.8, b=0.93, c=-0.169):
    """
    kdp = Specific Differential Phase (deg km^-1)
    zdr = Differential Reflectivity (dB)
    a, b, c = Adjustable function parameters
    """
    return a * kdp**b * 10.0**(c * zdr)


def calc_rain_z_zdr(dz, zdr, a=6.7e-3, b=0.927, c=-0.343):
    """
    dz = Reflectivity (dBZ)
    zdr = Differential Reflectivity (dB)
    a, b, c = Adjustable function parameters
    """
    return a * linearize(dz)**b * 10.0**(c * zdr)


def calc_rain_kdp(kdp, a=40.5, b=0.85):
    """
    kdp = Specific Differential Phase (deg km^-1)
    a, b = Adjustable coefficient, exponent
    """
    return a * kdp**b


def calc_blended_rain(
        dz=None, zdr=None, kdp=None, ice_flag=False, band='S',
        thresh_dz=38.0, thresh_zdr=0.5, thresh_kdp=0.3,
        thresh_frac_ice=0.1, thresh_nexrad=53.0,
        r_z_a=300.0, r_z_b=1.4, r_kdp_a=40.5, r_kdp_b=0.85,
        r_z_zdr_a=6.7e-3, r_z_zdr_b=0.927, r_z_zdr_c=-0.343,
        fit_a=None, fit_b=None, method='cr1995',
        r_kdp_zdr_a=90.8, r_kdp_zdr_b=0.93, r_kdp_zdr_c=-0.169):

    """
    This algorithm ingests polarimetric radar data and computes rain rate,
    using difference reflectivity as guide for determining presence of ice.

    Inputs:
    dz = Reflectivity
    zdr = Differential Reflectivity
    kdp = Specific Differential Phase
    thresh_dz = Threshold for dz to use certain rain algorithms
    thresh_zdr = Threshold for zdr to use certain rain algorithms
    thresh_kdp = Threshold for kdp to use certain rain algorithms
    thresh_nexrad = Threshold for dz to cap Z-R
    r_z_a, r_z_b => Z = a * R**b
    r_kdp_a, r_kdp_b => R = a * kdp**b
    r_kdp_zdr_a, r_kdp_zdr_b, r_kdp_zdr_c => R = a * kdp**b * 10.0**(c * zdr)
    r_z_zdr_a, r_z_zdr_b, r_z_zdr_c => R = a*linearize(z)**b * 10.0**(c*zdr)
    ice_flag = Set to True to return Zdp and Fraction of Ice
    fit_a, fit_b = Parameters for the rain line
    method = Reference used to set the rain line, if fit_a/b not provided

    Returns: rain rate, method, (Zdp, Fraction of Ice)
    method = 1: R(Kdp, Zdr)
    method = 2: R(Kdp)
    method = 3: R(Z, Zdr)
    method = 4: R(Z)
    method = 5: R(Zrain)

    See Bringi and Chandrasekar textbook for more information
    """
    # Initialize, check for all vars, check for scalars
    if dz is None or kdp is None or zdr is None:
        warnings.warn('No dz, zdr, or kdp provided, failing ...')
        return
    dz, zdr, kdp, len_flag = _check_for_array(dz, zdr, kdp)
    r_blended = np.zeros_like(dz)
    meth = np.int16(np.zeros_like(dz))

    # Rainfall estimates regardless of frequency
    r_dz_nexrad = calc_rain_nexrad(dz, thresh_nexrad, a=r_z_a, b=r_z_b)
    r_kdp = calc_rain_kdp(kdp, a=r_kdp_a, b=r_kdp_b)

    # Conditions used regardless of frequency
    cond_kdp = kdp >= thresh_kdp
    cond_dz = dz >= thresh_dz
    cond_dz_kdp = np.logical_and(cond_dz, cond_kdp)

    if band == 'S':
        # S-band only rain calculations
        r_kdp_zdr = calc_rain_kdp_zdr(kdp, zdr, a=r_kdp_zdr_a,
                                      b=r_kdp_zdr_b, c=r_kdp_zdr_c)
        r_dz_zdr = calc_rain_z_zdr(dz, zdr, a=r_z_zdr_a,
                                   b=r_z_zdr_b, c=r_z_zdr_c)

        # calculate rain line, zdp, and ice fraction
        zhor, zvert = calc_zh_zv(dz, zdr)
        zdp = calc_zdp(zhor, zvert)

        # calculate contribution to Zh from pure rain
        if fit_a is None:
            fit_a, fit_b = get_linear_fits(method=method)
        zrain = 10.0**((fit_a * zdp + fit_b)/10.0)
        dzrain = 10.0 * np.log10(zrain)
        fi = 1.0 - (zrain / zhor)
        r_dz_rainonly = calc_rain_nexrad(dzrain, thresh_nexrad,
                                         a=r_z_a, b=r_z_b)

        # Additional S-band conditions
        cond_ice = fi < thresh_frac_ice
        cond_zdr = zdr >= thresh_zdr

        # Set of method choices
        cond_meth_1or2 = np.logical_and(cond_ice, cond_dz_kdp)
        cond_meth_1 = np.logical_and(cond_meth_1or2, cond_zdr)
        cond_meth_2a = np.logical_and(cond_meth_1or2, ~cond_zdr)
        cond_meth_2b = np.logical_and(~cond_ice, cond_dz_kdp)
        cond_meth_2 = np.logical_or(cond_meth_2a, cond_meth_2b)
        cond_meth_3or4 = np.logical_and(cond_ice, ~cond_dz_kdp)
        cond_meth_3 = np.logical_and(cond_meth_3or4, cond_zdr)
        cond_meth_4 = np.logical_and(cond_meth_3or4, ~cond_zdr)
        cond_meth_5 = np.logical_and(~cond_ice, ~cond_dz_kdp)

        # Assign methods native only to S-band
        meth[cond_meth_1] = 1
        r_blended[cond_meth_1] = r_kdp_zdr[cond_meth_1]
        meth[cond_meth_3] = 3
        r_blended[cond_meth_3] = r_dz_zdr[cond_meth_3]
        meth[cond_meth_5] = 5
        r_blended[cond_meth_5] = r_dz_rainonly[cond_meth_5]

    else:
        # Logic for other bands; Zdr not used.
        # This will get spoofed by pure hail shaft (high Z, low Kdp)
        # Recommend csu_hidro_rain if hail is possible
        cond_meth_2 = cond_dz_kdp
        cond_meth_4 = ~cond_dz_kdp

    # Assign methods native to all bands
    meth[cond_meth_2] = 2
    r_blended[cond_meth_2] = r_kdp[cond_meth_2]
    meth[cond_meth_4] = 4
    r_blended[cond_meth_4] = r_dz_nexrad[cond_meth_4]

    # Return based on what the user provided and what they wanted
    if not ice_flag or band != 'S':
        if not len_flag:
            return r_blended[0], meth[0]
        else:
            return r_blended, meth
    else:
        if not len_flag:
            return r_blended[0], meth[0], zdp[0], fi[0]
        else:
            return r_blended, meth, zdp, fi


def csu_hidro_rain(
        dz=None, zdr=None, kdp=None, fhc=None, band='S',
        thresh_dz=38.0, thresh_zdr=0.5, thresh_kdp=0.3, thresh_nexrad=53.0,
        r_z_a=300.0, r_z_b=1.4, r_kdp_a=40.5, r_kdp_b=0.85,
        r_z_zdr_a=6.7e-3, r_z_zdr_b=0.927, r_z_zdr_c=-0.343,
        r_kdp_zdr_a=90.8, r_kdp_zdr_b=0.93, r_kdp_zdr_c=-0.169):
    """
    This algorithm ingests polarimetric radar data and computes rain rate,
    using hydrometeor ID as a guide for determining the presence of ice.

    Inputs:
    dz = Reflectivity
    zdr = Differential Reflectivity
    kdp = Specific Differential Phase
    fhc = Hydrometeor ID
    thresh_dz = Threshold for dz to use certain rain algorithms
    thresh_zdr = Threshold for zdr to use certain rain algorithms
    thresh_kdp = Threshold for kdp to use certain rain algorithms
    thresh_nexrad = Threshold for dz to cap Z-R
    r_z_a, r_z_b => Z = a * R**b
    r_kdp_a, r_kdp_b => R = a * kdp**b
    r_kdp_zdr_a, r_kdp_zdr_b, r_kdp_zdr_c => R = a * kdp**b * 10.0**(c * zdr)
    r_z_zdr_a, r_z_zdr_b, r_z_zdr_c => R = a*linearize(z)**b * 10.0**(c*zdr)

    Returns: rain rate, method
    method = 1: R(Kdp, Zdr)
    method = 2: R(Kdp)
    method = 3: R(Z, Zdr)
    method = 4: R(Z)

    See Bringi and Chandrasekar textbook for more information
    fixed 9/16/2015 - No rain was being calculated when low Z & high Kdp/Zdr
                      yet cond_rain == True
    """
    # Initialize, check for all necessary vars, and allow scalars
    if dz is None or kdp is None or zdr is None or fhc is None:
        warnings.warn('No dz, zdr, kdp, or fhc provided, failing ...')
        return
    dz, zdr, kdp, len_flag = _check_for_array(dz, zdr, kdp)
    crr_hidzk = np.zeros_like(dz)
    crr_meth = np.int16(np.zeros_like(dz))

    # Rainfall calculations regardless of frequency
    crr_z88t = calc_rain_nexrad(dz, a=r_z_a, b=r_z_b)
    crr_kd = calc_rain_kdp(kdp, a=r_kdp_a, b=r_kdp_b)

    # Conditions used regardless of frequency
    # fhc = 1, 2, 10 => drizzle, rain, big drops
    cond_rain_a = np.logical_or(fhc == 1, fhc == 2)
    cond_rain = np.logical_or(fhc == 10, cond_rain_a)
    cond_kdp = kdp >= thresh_kdp
    cond_dz = dz >= thresh_dz
    cond_dz_kdp = np.logical_and(cond_dz, cond_kdp)
    # fhc = 5, 7, 8, 9 => wet snow, graupel (LD & HD), hail
    cond_ice_a = np.logical_or(fhc == 5, fhc == 7)
    cond_ice_b = np.logical_or(fhc == 8, fhc == 9)
    cond_ice = np.logical_or(cond_ice_a, cond_ice_b)
    cond_ice_kdp = np.logical_and(cond_ice, cond_kdp)

    if band == 'S':
        # S-band only rain calcs
        crr_zhdr = calc_rain_z_zdr(dz, zdr, a=r_z_zdr_a,
                                   b=r_z_zdr_b, c=r_z_zdr_c)
        crr_kddr = calc_rain_kdp_zdr(kdp, zdr, a=r_kdp_zdr_a,
                                     b=r_kdp_zdr_b, c=r_kdp_zdr_c)

        # S-band only conditions
        cond_zdr = zdr >= thresh_zdr
        cond_dz_kdp_zdr = np.logical_and(cond_dz_kdp, cond_zdr)
        cond_dz_kdp_not_zdr = np.logical_and(cond_dz_kdp, ~cond_zdr)
        cond_not_kdp_zdr = np.logical_and(~cond_kdp, cond_zdr)
        cond_not_kdp_not_zdr = np.logical_and(~cond_kdp, ~cond_zdr)

        # Set methods
        cond_meth_1 = np.logical_and(cond_rain, cond_dz_kdp_zdr)
        cond_meth_2 = np.logical_and(cond_rain, cond_dz_kdp_not_zdr)
        cond_meth_3 = np.logical_and(cond_rain, cond_not_kdp_zdr)
        cond_meth_4a = np.logical_and(cond_rain, cond_not_kdp_not_zdr)
        cond_meth_4b = np.logical_and(cond_rain, ~cond_dz)
        cond_meth_4 = np.logical_or(cond_meth_4a, cond_meth_4b)

        # Assign rain rates based on methods
        crr_meth[cond_meth_1] = 1
        crr_hidzk[cond_meth_1] = crr_kddr[cond_meth_1]
        crr_meth[cond_meth_3] = 3
        crr_hidzk[cond_meth_3] = crr_zhdr[cond_meth_3]
        crr_meth[cond_ice_kdp] = 2
        crr_hidzk[cond_ice_kdp] = crr_kd[cond_ice_kdp]

    else:
        # Logic for other bands; Zdr not used.
        cond_meth_2 = cond_dz_kdp
        cond_meth_4 = np.logical_and(~cond_dz_kdp, cond_rain)
        # (high Z, low Kdp) or (low Z, high Kdp) and no HID rain, R = 0.0
        cond_ice_kdp = np.logical_and(~cond_dz_kdp, ~cond_rain)
        crr_meth[cond_ice_kdp] = 2
        crr_hidzk[cond_ice_kdp] = 0.0

    # Assign methods native to all bands
    crr_meth[cond_meth_2] = 2
    crr_hidzk[cond_meth_2] = crr_kd[cond_meth_2]
    crr_meth[cond_meth_4] = 4
    crr_hidzk[cond_meth_4] = crr_z88t[cond_meth_4]

    if not len_flag:
        return crr_hidzk[0], crr_meth[0]
    else:
        return crr_hidzk, crr_meth


def _check_for_array(dz, zdr, kdp):
    len_flag = hasattr(dz, '__len__')
    if not len_flag:
        dz = np.array([dz])
        kdp = np.array([kdp])
        zdr = np.array([zdr])
    return dz, zdr, kdp, len_flag
