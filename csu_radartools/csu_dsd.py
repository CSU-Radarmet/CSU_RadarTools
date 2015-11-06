"""
#############
csu_dsd sub-module of csu_radartools

Contacts
--------
Brenda Dolan (dolan@atmos.colostate.edu)
Timothy Lang (tjlangco@gmail.com)

References
----------
Bringi et al. (2004; JTECH) - Alternate S-band retrieval
Bringi et al. (2009; JTECH) - C-band retrieval
Bringi et al. (2013; AMS Radar Conf) - S-band retrieval

#############
"""

from __future__ import absolute_import
from __future__ import division
import numpy as np
from warnings import warn
from .csu_liquid_ice_mass import linearize
from .csu_blended_rain import _check_for_array

DEFAULT_MU = 3.0


def calc_dsd(dz=None, zdr=None, kdp=None, band='S', method='2013'):
    """This is the primary function to use, it calls everything else"""
    # Initialize, check for all vars, check for scalars
    if dz is None or kdp is None or zdr is None:
        warnings.warn('No dz, zdr, or kdp provided, failing ...')
        return
    dz, zdr, kdp, len_flag = _check_for_array(dz, zdr, kdp)
    # Get DSD based on wavelength
    if band == 'S':
        if method == '2013':
            d0, Nw, mu = calc_dsd_sband_bringi_2013(dz, zdr)
        else:
            d0, Nw, mu = calc_dsd_sband_bringi_2004(dz, zdr, kdp)
    elif band == 'C':
        d0, Nw, mu = calc_dsd_cband_bringi_2009(dz, zdr)
    else:
        warn('Unknown or unsupported wavelength band, try C or S')
        return
    if not len_flag:
        return d0[0], Nw[0], mu[0]
    else:
        return d0, Nw, mu


def power_law(var, a=1.0, b=1.0):
    """Generic power law equation"""
    return a * var**b

######################################################
# Bringi et al. 2009 retrievals
# Could make the following a class, but prefer to keep
# easy independent access to sub-functions
######################################################


def calc_dsd_cband_bringi_2009(dz, zdr):
    """
    Retrieves d0, Nw, and mu following the methodology of
    Bringi et al. (2009)
    Works for C-band
    """
    d0 = 0.0 * dz
    Nw = 0.0 * dz
    dz_lin = linearize(dz)
    meth1 = np.logical_and(zdr >= -0.5, zdr <= 1.25)
    d0[meth1] = d0_low_zdr(zdr[meth1])
    meth2 = zdr > 1.25
    d0[meth2] = d0_high_zdr(zdr[meth2])
    Nw = dz_lin / power_law(d0, 0.056, 7.319)
    mu = 0.0 * dz + DEFAULT_MU
    return d0, Nw, mu


def d0_low_zdr(zdr, a=0.0203, b=-0.1488, c=0.2209, d=0.5571, e=0.801):
    return power_law(zdr, a, 4) + power_law(zdr, b, 3) + \
        power_law(zdr, c, 2) + power_law(zdr, d) + e


def d0_high_zdr(zdr, a=0.0355, b=-0.3021, c=1.0556, d=0.6844):
    return power_law(zdr, a, 3) + power_law(zdr, b, 2) + power_law(zdr, c) + d

######################################################
# Bringi et al. 2013 retrievals
# Leverages Bringi et al. (2009) function design
######################################################


def calc_dsd_sband_bringi_2013(dz, zdr):
    """
    Retrieves d0, Nw, and mu following the methodology of
    Bringi et al. (2013)
    Works for S-band
    """
    d0 = 0.0 * dz
    Nw = 0.0 * dz
    dz_lin = linearize(dz)
    meth1 = zdr < 1
    d0[meth1] = d0_low_zdr(zdr[meth1], a=0.0424, b=-0.4571, c=0.6125, d=0.457,
                           e=0.8808)
    meth2 = zdr >= 1
    d0[meth2] = d0_high_zdr(zdr[meth2], a=0.0536, b=-0.1971,
                            c=0.6261, d=1.0815)
    Nw = 19.76 * dz_lin / power_law(d0, b=7.46)
    mu = 0.0 * dz + DEFAULT_MU
    return d0, Nw, mu

######################################################
# Bringi et al. 2004 retrievals
# Could make the following a class, but prefer to keep
# easy independent access to sub-functions
######################################################


def calc_dsd_sband_bringi_2004(dz, zdr, kdp):
    """
    Retrieves d0, Nw, and mu following the methodology of
    Bringi et al. (2004)
    Works for S-band
    """
    d0 = 0.0 * dz
    Nw = 0.0 * dz
    mu = 0.0 * dz
    zeta = linearize(zdr)
    dz_lin = linearize(dz)
    beta = calc_beta(dz_lin, kdp, zeta)
    # Break down DSD calculations by methods
    meth1a = np.logical_and(dz >= 35, zdr >= 0.2)
    meth1 = np.logical_and(meth1a, kdp >= 0.3)
    d0[meth1], Nw[meth1], mu[meth1] = dsd_sband_method1(
        dz_lin[meth1], kdp[meth1], zeta[meth1], beta[meth1])
    meth2 = np.logical_and(kdp < 0.3, zdr >= 0.5)
    d0[meth2], Nw[meth2], mu[meth2] = dsd_sband_method2(
        dz_lin[meth2], kdp[meth2], zdr[meth2])
    meth3a = np.logical_and(zdr >= 0.2, zdr < 0.5)
    meth3 = np.logical_and(meth3a, kdp < 0.3)
    d0[meth3], Nw[meth3], mu[meth3] = dsd_sband_method3(
        dz_lin[meth3], kdp[meth3], zdr[meth3], zeta[meth3])
    meth4a = np.logical_and(zdr >= -0.5, zdr < 0.2)
    meth4 = np.logical_and(meth4a, kdp < 0.3)
    d0[meth4], Nw[meth4], mu[meth4] = dsd_sband_method4(dz_lin[meth4],
                                                        zeta[meth4])
    return d0, Nw, mu


def calc_beta(dz, kdp, zeta, a=2.08, b=-0.365, c=0.38, d=0.965):
    """Bringi et al. (2004) Eq. A.1"""
    return a * (dz**b) * (kdp**c) * (zeta**d)


def calc_d0_sband_method1(dz, zeta, beta, a=0.56, b=0.064, c=0.024, d=-1.42):
    """Bringi et al. (2004) Eqs. A.2-A.5 """
    c1 = power_law(beta, c, d)
    return d0_from_dz_zeta(dz, zeta, a=a, b=b, c=c1)


def calc_logNw_sband_method1(dz, zeta, beta, a=3.29, b=0.058, c=-0.023,
                             d=-1.389):
    """Bringi et al. (2004) Eqs. A.6-A.9 """
    return a * (dz**b) * (zeta**(power_law(beta, c, d)))


def calc_mu_sband_method1(d0, zeta, beta, a1=203.0, a2=1.89,
                          b1=2.23, b2=0.0388,
                          c1=3.16, c2=-0.0463, d1=0.374, d2=-0.355):
    """Bringi et al. (2004) Eqs. A.10-A.14 """
    a = power_law(beta, a1, a2)
    b = power_law(beta, b1, b2)
    c = power_law(beta, c1, c2)
    d = power_law(beta, d1, d2)
    return a * ((d0**b) / (zeta-1.0)) - power_law(zeta, c, d)


def calc_d0_sband_method3(dz, zdr, zeta, a=1.81, b=0.486, c=0.6096,
                          d=0.0516, e=3.111):
    """Bringi et al. (2004) Eq. A.18"""
    part1 = ((zdr-0.2)/3.0) * power_law(zdr, a, b)
    part2 = ((0.5-zdr)/3.0) * d0_from_dz_zeta(dz, zeta)
    return part1 + part2


def calc_Nw_mult_methods(dz, d0):
    """Bringi et al. (2004) Eqs. A.17, A.19, A.21 """
    return 21.0 * dz / power_law(d0, b=7.3529)


def d0_from_dz_zeta(dz, zeta, a=0.6096, b=0.0516, c=3.111):
    """Bringi et al. (2004) Eq. A.20"""
    return a * dz**b * zeta**c


def dsd_sband_method1(dz, kdp, zeta, beta):
    """Bringi et al. (2004) Eqs. A.1-A.14 """
    d0 = calc_d0_sband_method1(dz, zeta, beta)
    logNw = calc_logNw_sband_method1(dz, zeta, beta)
    mu = calc_mu_sband_method1(d0, zeta, beta)
    return d0, 10.0**logNw, mu


def dsd_sband_method2(dz, kdp, zdr):
    """Bringi et al. (2004) Eqs. A.15-A.17 """
    d0 = power_law(zdr, 1.81, 0.486)
    Nw = calc_Nw_mult_methods(dz, d0)
    mu = 0.0 * dz + DEFAULT_MU
    return d0, Nw, mu


def dsd_sband_method3(dz, kdp, zdr, zeta):
    """Bringi et al. (2004) Eqs. A.18-A.19 """
    d0 = calc_d0_sband_method3(dz, zdr, zeta)
    Nw = calc_Nw_mult_methods(dz, d0)
    mu = 0.0 * dz + DEFAULT_MU
    return d0, Nw, mu


def dsd_sband_method4(dz, zeta):
    """Bringi et al. (2004) Eqs. A.20-A.21 """
    d0 = d0_from_dz_zeta(dz, zeta)
    Nw = calc_Nw_mult_methods(dz, d0)
    mu = 0.0 * dz + DEFAULT_MU
    return d0, Nw, mu

######################################################
