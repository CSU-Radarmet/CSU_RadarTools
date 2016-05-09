# -*- coding: utf-8 -*-
"""
attenuation sub-module of csu_radartools

Functions to calculate coefficients used with attenuation calculations.

References
----------
Battan (1973), Radar Observations of the Atmosphere.
Doviak and Zrnic (1993), Doppler Radar and Weather Observations.
"""
import numpy as np


def abs_coeff(D, lam, m):
    """
    Absorption coefficient of a spherical particle. Unitless.

    Doviak and Zrnic (1993), Eqn 3.14a or Battan (1973), Eqn 6.6

    Parameters
    ----------
    D : float or array
        Particle diameter [m]
    lam : float
        Radar wavelength [m]
    m : float
        Complex refractive index [unitless], in the form 7.14 - 2.89j

    Notes
    -----
    An example from Battan (1973) is for water at 0C m=7.14-2.89j for a
       wavelength of 3.21 cm and for ice m=1.78-0.0024j for
       wavelength range from 1-10 cm.
    See Battan (1973) Ch.4 , Tables 4.1 and 4.2 for values from
       Gunn and East (1954).
    Also see Doviak and Zrnic (1993), Fig. 3.3 caption.
    """
    Km = (m**2 - 1) / (m**2 + 2)
    Qa = (np.pi**2 * np.asarray(D)**3 / lam) * np.imag(-1 * Km)
    return Qa


def scat_coeff(D, lam, m):
    """
    Scattering coefficient of a spherical particle. Unitless.

    Doviak and Zrnic (1993), Eqn 3.14b or Battan (1973), Eqn 6.5

    Parameters
    ----------
    D : float or array
        Particle diameter [m]
    lam : float
        Radar wavelength [m]
    m : float
        Complex refractive index [unitless]

    Notes
    -----
    An example from Battan (1973) is for water at 0C m=7.14-2.89j for a
       wavelength of 3.21 cm and for ice m=1.78-0.0024j for
       wavelength range from 1-10 cm.
    See Battan (1973) Ch.4 , Tables 4.1 and 4.2 for values from
       Gunn and East (1954).
    Also see Doviak and Zrnic (1993), Fig. 3.3 caption.
    """

    Km = (m**2 - 1) / (m**2 + 2)
    Qs = (2 * np.pi**5 * np.asarray(D)**6 / (3 * lam**4) * (np.absolute(Km))**2)
    return Qs


def ext_coeff(D, lam, m):
    """
    Extinction coefficient of a spherical particle. Unitless.

    Doviak and Zrnic (1993), Eqn 3.14b or Battan (1973), Eqn 6.5

    Parameters
    ----------
    D : float or array
        Particle diameter [m]
    lam : float
        Radar wavelength [m]
    m : float
        Complex refractive index [unitless]

    USAGE::
    -----
    Qe = ext_coeff(D,lam,m)

    NOTES::
    -----
    An example from Battan (1973) is for water at 0C m=7.14-2.89j for a
       wavelength of 3.21 cm and for ice m=1.78-0.0024j for
       wavelength range from 1-10 cm.
    See Battan (1973) Ch.4 , Tables 4.1 and 4.2 for values from
       Gunn and East (1954).
    Also see Doviak and Zrnic (1993), Fig. 3.3 caption.
    """

    Qa = abs_coeff(np.asarray(D), lam, m)
    Qs = scat_coeff(np.asarray(D), lam, m)
    Qe = Qa + Qs

    return Qe
