# -*- coding: utf-8 -*-
"""
csu_blended_rain_tropical.py

# Brody Fuchs, CSU, Oct 2014
# brfuchs@atmos.colostate.edu

# python remake of Brenda's code to calculate water and ice mass

Amendments by
Timothy Lang (tjlangco@gmail.com)
2/20/2015

Designed around tropical, oceanic equations found in Thompson et al. 2016, and
designed to work with the Powell and Houze 2015 rain typing (or simple
convective /stratiform)
Brenda Dolan (bdolan@atmos.colostate.edu)
8/2016

Added the dBZ + Kdp threshold back in to eliminate small reflectivities with non-zero Kdp
blowing up the rain rates.
11/2016
"""

from __future__ import absolute_import
from __future__ import division
import numpy as np
import warnings
from .csu_liquid_ice_mass import linearize, get_linear_fits
from .common import (
    _check_for_array, calc_rain_zr, calc_rain_nexrad, calc_rain_kdp_zdr,
    calc_rain_z_zdr, calc_rain_kdp)


def calc_blended_rain_tropical(
        dz=None, zdr=None, kdp=None, cs=None, ice_flag=False,
        fhc=None,
        predef='True', band='S',
        thresh_dz=38.0, thresh_zdr=0.25, thresh_kdp=0.3,
        thresh_frac_ice=0.1, thresh_nexrad=53.0,
        r_z_a=216., r_z_b=1.39,
        r_z_a_c=126., r_z_b_c=1.39, r_z_a_s=291., r_z_b_s=1.55,
        r_kdp_a=59.5202, r_kdp_b=0.7451,
        r_z_zdr_a=0.0085, r_z_zdr_b=0.9237, r_z_zdr_c=-5.2389/10.,
        fit_a=None, fit_b=None, method='cr1995',
        r_kdp_zdr_a=96.5726, r_kdp_zdr_b=0.9315, r_kdp_zdr_c=-2.1140/10.):

    """
    This algorithm ingests polarimetric radar data and computes rain rate,
    based on Thompson et al. 2016. Since ice is not expected at the surface in
    the tropics, the algorithm does not look for ice contamination.

    Inputs:
    dz = Reflectivity
    zdr = Differential Reflectivity
    kdp = Specific Differential Phase
    cs = Convective / Stratiform map (
        3=mixed / uncertain,
        2=convective,
        1= stratiform,
        0=unknown)
    thresh_zdr = Threshold for zdr to use certain rain algorithms
    thresh_kdp = Threshold for kdp to use certain rain algorithms

    Use the predef='True' if you want to use the coefficients defined
    in Thompson 2016. Otherwise specifiy the coefficients for each
    relationship.

    Returns: rain rate, method, (Zdp, Fraction of Ice)
    method = 1: R(Kdp, Zdr)
    method = 2: R(Kdp)
    method = 3: R(Z, Zdr)
    method = 4: R(Z_all)
    method = 5: R(Z_c)
    method = 6: R(Z_s)

    Implemented by B. Dolan July 2016

    See Bringi and Chandrasekar textbook for more information
    See also Thompson et al. 2016 for more details of the Tropical Blended
    algorithm
    """
    bnd = band
    # Initialize, check for all vars, check for scalars
    if dz is None or kdp is None or zdr is None:
        warnings.warn('No dz, zdr, or kdp provided, failing ...')
        return
    else:
        len_flag = hasattr(dz, '__len__')

        if not len_flag:
            dz = np.array([dz])
            kdp = np.array([kdp])
            zdr = np.array([zdr])

    """
    NOTE: The T15 equations are defined as:
        R = a*Kdp*zeta_dr**c and R=a*Zh*zeta_dr**c
        where as the algorithm expects it in the form
        R = a*Kdp*10.**(Zdr*c)  R=a*Zh*10.**(Zdr*c)

        These just differ in the by 1/10. in the exponent of Zdr, so be sure to
        divide by 10. when sending to the algorithm.
    """
    if predef == 'True':
        """
        NOTE: The coefficients are defined as:
        Z=aR**b
        """
        # Convective R-Z
        r_z_a_c = 126.
        r_z_b_c = 1.39
        # Stratiform R-Z
        r_z_a_s = 291.
        r_z_b_s = 1.55
        # All R-Z
        r_z_a = 216
        r_z_b = 1.39

        if band == 'S':
            # R = a*Kdp*zeta**c
            r_kdp_zdr_a = 96.5726,
            r_kdp_zdr_b = 0.9315,
            r_kdp_zdr_c = -2.1140/10.,

            # R = a*Zh*zeta**c
            r_z_zdr_a = 0.0085
            r_z_zdr_b = 0.9237
            r_z_zdr_c = -5.2389/10.

            # R = a*Kdp**b
            r_kdp_a = 59.5202
            r_kdp_b = 0.7451

        if band == 'C':
            # R = a*Kdp*zeta**c
            r_kdp_zdr_a = 45.6976
            r_kdp_zdr_b = 0.8763
            r_kdp_zdr_c = -1.6718/10.

            # R = a*Zh*zeta**c
            r_z_zdr_a = 0.0086
            r_z_zdr_b = 0.9088
            r_z_zdr_c = -4.2059/10.

            # R = a*Kdp**b
            r_kdp_a = 34.5703
            r_kdp_b = 0.7331

        if band == 'X':
            # R = a*Kdp*zeta**c
            r_kdp_zdr_a = 28.1289
            r_kdp_zdr_b = 0.9194
            r_kdp_zdr_c = -1.6876/10.

            # R = a*Zh*zeta**c
            r_z_zdr_a = 0.0085
            r_z_zdr_b = 0.9294
            r_z_zdr_c = -4.4580/10.

            # R = a*Kdp**b
            r_kdp_a = 21.9729
            r_kdp_b = 0.7221

    if fhc is None:
        warnings.warn('No FHC ... Rain may be calculated above melting layer')

    r_blended = np.zeros_like(dz)
    # r_blended[...] = -9999.
    meth = np.int16(np.zeros_like(dz)) - 1

    # convective ZR estimate
    r_dz_conv = calc_rain_zr(dz, a=r_z_a_c, b=r_z_b_c)
    # stratiform ZR estimate
    r_dz_strat = calc_rain_zr(dz, a=r_z_a_s, b=r_z_b_s)
    # ALL ZR estimate
    r_dz_all = calc_rain_zr(dz, a=r_z_a, b=r_z_b)

    # Polarimetric estimates
    r_kdp_zdr = calc_rain_kdp_zdr(kdp, zdr, a=r_kdp_zdr_a,
                                  b=r_kdp_zdr_b, c=r_kdp_zdr_c)
    r_dz_zdr = calc_rain_z_zdr(dz, zdr, a=r_z_zdr_a,
                               b=r_z_zdr_b, c=r_z_zdr_c)
    r_kdp = calc_rain_kdp(kdp, a=r_kdp_a, b=r_kdp_b)

    # Conditions
    cond_kdp = kdp >= thresh_kdp
    cond_dbz = dz >= thresh_dz
    cond_kdpz = np.logical_and(cond_kdp,cond_dbz)
    cond_zdr = zdr >= thresh_zdr
    # Set of method choices
    cond_meth_1 = np.logical_and(cond_kdpz, cond_zdr)
    cond_meth_2 = np.logical_and(cond_kdpz, ~cond_zdr)
    cond_meth_3 = np.logical_and(cond_zdr, ~cond_kdpz)
    cond_meth_4 = np.logical_and(~cond_kdpz, ~cond_zdr)

    if cs is None:
        meth[cond_meth_4] = 4
        r_blended[cond_meth_4] = r_dz_all[cond_meth_4]
    else:
        cond_meth_4c = np.logical_and(cond_meth_4, cs == 2)
        cond_meth_4s = np.logical_and(cond_meth_4, cs == 1)
        cond_meth_4a = np.logical_and(cond_meth_4, cs == 3)
        meth[cond_meth_4c] = 5
        meth[cond_meth_4s] = 6
        meth[cond_meth_4a] = 4
        r_blended[cond_meth_4c] = r_dz_conv[cond_meth_4c]
        r_blended[cond_meth_4s] = r_dz_strat[cond_meth_4s]
        r_blended[cond_meth_4a] = r_dz_all[cond_meth_4a]

    # Assign methods
    meth[cond_meth_1] = 1
    meth[cond_meth_2] = 2
    meth[cond_meth_3] = 3
    # Assign rain rates based on methods
    r_blended[cond_meth_1] = r_kdp_zdr[cond_meth_1]
    r_blended[cond_meth_2] = r_kdp[cond_meth_2]
    r_blended[cond_meth_3] = r_dz_zdr[cond_meth_3]

    if fhc is not None:
        cond_ice = np.logical_and(fhc > 2, fhc < 10)
        r_blended[cond_ice] = 0.0
        meth[cond_ice] = -1
        cond_hail = np.logical_and(fhc == 9, cond_kdpz)
        r_blended[cond_hail] = r_kdp[cond_hail]
        meth[cond_ice] = 2

    r_blended[dz < -10] = 0.0
    meth[dz < -10] = -1

    # Return based on what the user provided and what they wanted
    return r_blended, meth
