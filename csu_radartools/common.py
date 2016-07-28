# -*- coding: utf-8 -*-
"""
common sub-module of csu_radartools

Contains commonly used functions.

Nick Guy 10 May 2016
Timothy Lang Rev1 (7/28/2016) - Moved functions from csu_blended_rain here.
"""
import numpy as np

##############
##  Arrays  ##
##############


def _check_for_array(dz, zdr, kdp):
    len_flag = hasattr(dz, '__len__')
    if not len_flag:
        dz = np.array([dz])
        kdp = np.array([kdp])
        zdr = np.array([zdr])
    return dz, zdr, kdp, len_flag


################
##  Rainfall  ##
################


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

########################
##  Unit Conversions  ##
########################


def dbz2z(dbz):
    """
    Convert from log [dBZ] to linear Z [mm^6 m^−3] units.

    Parameters
    ----------
    dbz : float or array
        logarithmic reflectivity value
    """
    return 10.**(np.asarray(dbz)/10.)


def linearize(dbz):
    return dbz2z(dbz)


def z2dbz(zlin):
    """
    Convert from linear Z [mm^6 m^−3] to log [dBZ] units.

    Parameters
    ----------
    zlin : float or array
        linear reflectivity units
    """
    return 10. * np.log10(np.asarray(zlin))


def si2kmh(si):
    """
    Convert from SI [m/s] wind units to km/h.

    Parameters
    ----------
    si : float or array
        Wind in SI units (m/s)
    """
    return np.asarray(si) * 3600. / 1000.


def si2mph(si):
    """
    Convert from SI wind units to miles/h [mph].

    Parameters
    ----------
    si: float or array
        Wind in SI units (m/s)
    """
    return np.asarray(si) * 0.62137 / 1000. * 3600.


def si2kts(si):
    """
    Convert from SI wind units to knots [kt].

    Parameters
    ----------
    si: float or array
        Wind in SI units (m/s)
    """
    return np.asarray(si) * 0.51


def kmh2si(kmh):
    """
    Convert from km/h to SI wind units [m/s].

    Parameters
    ----------
    kmh: float or array
        Wind in km/hr
    """
    return np.asarray(kmh) * 1000. / 3600.


def mph2si(mph):
    """
    Convert from miles/h to SI wind units [m/s].

    Parameters
    ----------
    mph: float or array
        Wind in miles per hour
    """
    return np.asarray(mph) * 1000. / (0.62137 * 3600.)


def kts2si(kts):
    """
    Convert from knots to SI wind units [m/s].

    Parameters
    ----------
    kts: float or array
        Wind in knots
    """
    return np.asarray(kts) / 0.51
