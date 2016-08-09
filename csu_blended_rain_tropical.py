"""
csu_blended_rain.py

# Brody Fuchs, CSU, Oct 2014
# brfuchs@atmos.colostate.edu

# python remake of Brenda's code to calculate water and ice mass

Amendments by
Timothy Lang (tjlangco@gmail.com)
2/20/2015

"""

import numpy as np
import warnings
from csu_liquid_ice_mass import linearize, get_linear_fits


def calc_rain_zr(dz, a=216, b=1.39, band='S'):
    """
    dz = Reflectivity (dBZ), returns rainfall rate (mm h^-1)
    Form Z = a * R**b
    """
    return (linearize(dz) / a)**(1.0 / b)


def calc_trop_rain_kdp_zdr(
    kdp,
    zdr,
    a=96.5726,
    b=0.315,
    c=-2.1140,
    band='S',
        predef='True'):
    """
    kdp = Specific Differential Phase (deg km^-1)
    zdr = Differential Reflectivity (dB)
    a, b, c = Adjustable function parameters
    """

    if predef == 'True':
        if band == 'S':
            a = 96.5726
            b = 0.315
            c = -2.1140

        if band == 'C':
            a = 45.6976
            b = 0.8763
            c = -1.6718

        if band == 'X':
            a = 28.1289
            b = 0.9194
            c = -1.6876

    return a * kdp**b * linearize(zdr)**c


def calc_trop_rain_z_zdr(
    dz,
    zdr,
    a=0.0085,
    b=0.9237,
    c=-5.2389,
    band='S',
        predef='True'):
    """
    dz = Reflectivity (dBZ)
    zdr = Differential Reflectivity (dB)
    a, b, c = Adjustable function parameters
    """
    if predef == 'True':
        if band == 'S':
            a = 0.0085
            b = 0.9237
            c = -5.2389

        if band == 'C':
            a = 0.0086
            b = 0.9088
            c = -4.2059

        if band == 'X':
            a = 0.0085
            b = 0.9294
            c = -4.4580

    return a * linearize(dz)**b * linearize(zdr)**c


def calc_trop_rain_kdp(kdp, a=59.5202, b=0.7451, band='S', predef='True'):
    """
    kdp = Specific Differential Phase (deg km^-1)
    a, b = Adjustable coefficient, exponent

    The coefficients below are from Thompson et al. 2016. They are for the
    convective regime since Kdp has to be significant enough to use, which is
    essentially for convective points.
    """
    if predef == 'True':
        if band == 'S':
            a = 59.5202
            b = 0.7451
        if band == 'C':
            a = 34.5703
            b = 0.7331
        if band == 'X':
            a = 21.9729
            b = 0.7221

    return a * kdp**b


def calc_blended_rain_tropical(
        dz=None,
        zdr=None,
        kdp=None,
        cs=None,
        thresh_zdr=0.25,
        thresh_kdp=0.3,
        fhc=None,
        band='S'):
    """
    This algorithm ingests polarimetric radar data and computes rain rate,
    based on Thompson et al. 2016. Since ice is not expected at the surface in
    the tropics, the algorithm does not look for ice contamination.

    Inputs:
    dz = Reflectivity
    zdr = Differential Reflectivity
    kdp = Specific Differential Phase
    cs = Convective / Stratiform map (2=convective, 1= unknown, 0=stratiform)
    thresh_zdr = Threshold for zdr to use certain rain algorithms
    thresh_kdp = Threshold for kdp to use certain rain algorithms

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

    if fhc is None:
        warnings.warn('No FHC ... Rain may be calculated above melting layer')

        r_blended = np.zeros_like(dz)
        r_blended[...] = -9999.
        meth = np.int16(np.zeros_like(dz))
        meth[...] = -1

        # convective ZR estimate
        r_dz_conv = calc_rain_zr(dz, a=126, b=1.39)
        # stratiform ZR estimate
        r_dz_strat = calc_rain_zr(dz, a=291, b=1.55)
        # ALL ZR estimate
        r_dz_all = calc_rain_zr(dz, a=216, b=1.39)

        # Polarimetric estimates
        r_dz_zdr = calc_trop_rain_z_zdr(dz, zdr, band=bnd, predef='True')
        r_kdp = calc_trop_rain_kdp(kdp, band=bnd, predef='True')
        r_kdp_zdr = calc_trop_rain_kdp_zdr(kdp, zdr, band=bnd, predef='True')

        # Conditions
        cond_kdp = kdp >= thresh_kdp
        cond_zdr = zdr >= thresh_zdr
        # Set of method choices
        cond_meth_1 = np.logical_and(cond_kdp, cond_zdr)
        cond_meth_2 = np.logical_and(cond_kdp, ~cond_zdr)
        cond_meth_3 = np.logical_and(cond_zdr, ~cond_kdp)
        cond_meth_4 = np.logical_and(~cond_kdp, ~cond_zdr)

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
            cond_ice = np.loglcal_or(fhc > 2, fhc < 10)
            r_blended[cond_ice] = -9999.
            meth[cond_ice] = -1

        r_blended[dz < -10] = -9999.
        meth[dz < -10] = -1

        # Return based on what the user provided and what they wanted
        return r_blended, meth


def _check_for_array(dz, zdr, kdp):
    len_flag = hasattr(dz, '__len__')
    if not len_flag:
        dz = np.array([dz])
        kdp = np.array([kdp])
        zdr = np.array([zdr])
    return dz, zdr, kdp, len_flag
