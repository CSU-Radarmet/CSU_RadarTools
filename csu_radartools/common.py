# -*- coding: utf-8 -*-
"""
common sub-module of csu_radartools

Contains commonly used functions.

Nick Guy 10 May 2016
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
