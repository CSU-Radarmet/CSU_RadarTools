# -*- coding: utf-8 -*-
"""
fundamentals sub-module of csu_radartools

Functions to calculate:
    radar system characteristics (System)
    Doppler radar system characteristics (Doppler)
    radar geometry characteristics (Geometry)
    radar-derived variables (Variables)
    coefficients used with attenuation calculations (Attenuation)

References
----------
Battan (1973), Radar Observations of the Atmosphere.
Doviak and Zrnic (1993), Doppler Radar and Weather Observations.
Rinehart (1997), Radar for Meteorologists.
Jorgensen (1983; JCAM), Feasibility Test of an Airborne Pulse-Doppler
    Meteorological Radar.
Bech et al. (2003; JAOT), The Sensitivity of Single Polarization Weather Radar
    Beam Blockage Correction to Variability in the Vertical Refractivity
    Gradient.
Aydin et al. (1986; JCAM), Remote Sensing of Hail with Dual Linear
    Polarization Radar

Nick Guy - 10 May 2016, Ported from PyRadarMet
"""
import numpy as np

speed_of_light = 3e8  # Speed of light [m/s]
k_boltz = 1.381e-23  # Boltzmann's constant [ m^2 kg s^-2 K^-1]
earth_radius = 6371000.  # Earth's average radius [m] assuming sphericity
r43 = earth_radius * 4./3.
# Earth radius according to International Union of Geodesy and Geophysics
# r43 is 4/3 Approximation effective radius for standard atmosphere [m]

############
#  System  #
############


def antenna_gain(p1, p2):
    """
    Antenna gain [dB] via power ratio.

    From Rinehart (1997), Eqn 2.1

    Parameters
    ----------
    p1 : float or array
        Power on the beam axis [W]
    p2 : float or array
        Power from an isotropic antenna [W]

    Notes
    -----
    Ensure that both powers have the same units!
    If arrays are used, either for one okay, or both must be the same length.
    """
    return 10. * np.log10(np.asarray(p1) / np.asarray(p2))


def frequency(wavelength):
    """
    Frequency [Hz] given wavelength.

    Parameters
    ----------
    wavelength : float or array
        Wavelength [m]
    """
    return speed_of_light / np.asarray(wavelength)


def wavelength(freq):
    """
    Wavelength [m] given frequency.

    Parameters
    ----------
    freq : float or array
        Frequency [Hz]
    """
    return speed_of_light / np.asarray(freq)


def pulse_length(pdur):
    """
    Pulse length [m] given pulse duration.

    Parameters
    ----------
    pDur : float or array
        Pulse duration [s]

    Notes
    -----
    This equation is only interested in return pulses to radar, leading
    to the factor of 1/2.
    """
    return speed_of_light * np.asarray(pdur) / 2.


def pulse_duration(tau):
    """
    Pulse duration [s] from pulse length.

    Parameters
    ----------
    tau : float or array
        Pulse length [m]
    """
    return 2 * np.asarray(tau) / speed_of_light


def radar_const(power_t, gain, tau, wavelength, bw_h, bw_v, aloss, rloss):
    """
    Radar constant. Unitless.

    From CSU Radar Meteorology notes, AT 741

    Parameters
    ----------
    power_t : float
        Transmitted power [W]
    gain : float
        Antenna Gain [dB]
    tau : float
        Pulse Width [s]
    wavelength : float
        Radar wavelength [m]
    bw_h : float
        Horizontalntenna beamwidth [degrees]
    bw_v : float
        Vertical antenna beamwidth [degrees]
    aloss : float
        Antenna/waveguide/coupler loss [dB]
    rloss : float
        Receiver loss [dB]
    """
    # Convert from dB to linear units
    alosslin = 10**(aloss / 10.)
    rlosslin = 10**(rloss / 10.)
    gainlin = 10**(gain / 10.)

    # Convert beamwidth to radians
    bw_hr = np.deg2rad(bw_h)
    bw_vr = np.deg2rad(bw_v)

    # Calculate the numerator
    numer = (np.pi**3 * speed_of_light * power_t * gainlin**2 * tau *
             bw_hr * bw_vr * alosslin * rlosslin)

    # Calculate the denominator
    denom = 1024. * np.log(2) * wavelength**2
    return numer/denom


def ant_eff_area(gain, wavelength):
    """
    Antenna effective area. [m^-2]

    From Rinehart (1997), Eqn 4.5

    Parameters
    ----------
    gain : float or array
        Antenna Gain [dB]
    wavelength : float
        Radar wavelength [m]
    """
    # Convert from dB to linear units
    gainlin = 10**(np.asarray(gain) / 10.)
    return gainlin * wavelength**2 / (4 * np.pi)


def power_target(power_t, gain, areat, range):
    """
    Power [W] intercepted by target.

    From Rinehart (1997), Eqn 4.3

    Parameters
    ----------
    power_t : float
        Transmitted power [W]
    gain : float
        Antenna gain [dB]
    areat : float
        Area of target [m^2]
    range : float or array
        Distance to sample volume from radar [m]
    """
    # Convert from dB to linear units
    gainlin = 10**(gain / 10.)
    return (power_t * gainlin * areat) / (4 * np.pi * np.asarray(range)**2)


def xsec_bscatter_sphere(diam, wavelength, dielectric=0.93):
    """
    Backscatter cross-sectional area [m^-2] of a sphere using
    the Rayleigh approximation.

    From Rinehart (1997), Eqn 4.9 and 5.7

    Parameters
    ----------
    diam : float or array
        Diamter of target [m]
    wavelength : float
        Radar wavelength [m]
    dielectric : float
        Dielectric factor [unitless]

    Notes
    -----
    The Rayleigh approximation is good when the diamter of a spherical particle
    is much smaller than the wavelength of the radar (D/wavelength= 1/16).
    This condition leads to the relationship that the area is proportional
    to the sixth power of the diameter.

    The default is for a dielectric factor value for water.  This can be
    changed by the user, e.g. K=0.208 for particle sizes of equivalent melted
    diameters or K=0.176 for particle sizes of equivalent ice spheres.
    """
    return (np.pi**5 * dielectric**2 * np.asarray(diam)**6) / wavelength**4


def norm_xsec_bscatter_sphere(diam, wavelength, dielectric=0.93):
    """
    Normalized Backscatter cross-sectional area [m^2] of a sphere using
    the Rayleigh approximation.

    From Rinehart (1997), Eqn 4.9 and 5.7 and Battan Ch. 4.5

    Parameters
    ----------
    diam : float or array
        Diamter of targer [m]
    wavelength : float
        Radar wavelength [m]
    dielectric : float
        Dielectric factor [unitless]

    Notes
    -----
    The Rayleigh approximation is good when the diamter of a spherical particle
    is much smaller than the wavelength of the radar (D/wavelength= 1/16).
    This condition leads to the relationship that the area is proportional
    to the sixth power of the diameter.

    The default is for a dielectric factor value for water.  This can be
    changed by the user, e.g. K=0.208 for particle sizes of equivalent melted
    diameters or K=0.176 for particle sizes of equivalent ice spheres.
    """

    # Calculate the cross-sectional backscatter area
    sig = xsec_bscatter_sphere(np.asarray(diam), wavelength, dielectric)
    return sig / (np.pi * (np.asarray(diam)/2.)**2)


def size_param(diam, wavelength):
    """
    Size parameter calculation. Unitless.

    From Rinehart (1997), Eqn 4.9 and 5.7 and Battan Ch. 4.5

    Parameters
    ----------
    diam : float or float
        Diamter of target [m]
    wavelength : float
        Radar wavelength [m]

    Notes
    -----
    The size paramter can be used along with the backscattering cross-section
    to distinguish ice and water dielectric characteristics.

    For example:
    Alpha < 2 the backscattering cross-section of ice is smaller than water,
    Alpha > 2 the opposite is true due to the fact that absorption in water
    exceeds that in ice.
    """
    return 2 * np.pi * np.asarray(diam)/2. / wavelength


def power_return_target(power_t, gain, wavelength, sig, range):
    """
    Power [W] returned y target located at the center of the antenna
    beam pattern.

    From Rinehart (1997), Eqn 4.7

    Parameters
    ----------
    power_t : float
        Transmitted power [W]
    gain : float
        Antenna gain [dB]
    wavelength : float
        Radar wavelength [m]
    sig : float
        Backscattering cross-sectional area of target [m^2]
    range : float or array
        Distance to sample volume from radar [m]
    """
    # Convert from dB to linear units
    gainlin = 10**(gain/10.)
    return ((power_t * gainlin**2 * wavelength**2 * sig) /
            (64 * np.pi**3 * np.asarray(range)**4))


def thermal_noise(bandwidth, units, noise_temp=290.):
    """
    Thermal noise power [W or 'dBm'].

    From CSU Radar Meteorology notes, AT741

    Parameters
    ----------
    bandwidth : float or array
        Receiver bandwidth [Hz]
    units : float
        String of nits desired, can be 'W' or 'dBm'
    noise_temp : float
        Reciever noise temperature [K]

    Notes
    -----
    Reciever noise temp set to conventional 290K by default
    """
    # Calculate the noise, convert if requested
    noise = k_boltz * noise_temp * np.asarray(bandwidth)

    if units.upper() == 'W':
        noiset = noise
    elif units.upper() == 'DBM':
        noiset = 10. * np.log10(noise/10**-3)
    else:
        print("Units must be in 'W' or 'dBm'")
        noiset = np.nan
    return noiset

#############
#  Doppler  #
#############


def fmax(prf):
    """Maximum frequency [Hz] given PRF.

    From Rinehart (1997), Eqn 6.8

    Parameters
    ----------
    prf : float or array
        Pulse repetition frequency [Hz]
    """
    return np.asarray(prf) / 2.


def vmax(prf, wavelength):
    """Nyquist velocity, or maximum unambiguous Doppler velocity (+ or -) [m/s].

    From Rinehart (1997), Eqn 6.7

    Parameters
    ----------
    prf : float or array
        Radar pulse repetition frequency [Hz]
    wavelength : float or array
        Radar wavelength [m]
    """
    return np.asarray(prf) * wavelength / 4.


def nyquist(prf, wavelength):
    """Wrapper function for vmax."""
    return vmax(prf, wavelength)


def rmax(prf):
    """Maximum unamiguous range [m].

    From Rinehart (1997), Eqn 6.11

    Parameters
    ----------
    prf : float or array
        Pulse repetition frequency [Hz]
    """
    return speed_of_light / (2. * np.asarray(prf))


def unambiguous_range(prf):
    """Wrapper for rmax."""
    return rmax(prf)


def doppler_dilemma(varin, wavelength):
    """
    The "Doppler dilemma" is the fact that both the Nyquist velocity and
    unambiguous maximum range of the radar are based upon the PRF of
    the system.

    However, they are inversely proportional, meaning that increasing one
    requires a decrease in the other.  A trade-off inherent in Doppler radar
    systems.  This relationship allows a solution for one variable given the
    value of the other.

    From Rinehart (1997), Eqn 6.12

    Parameters
    ----------
    varin : float or array
        Nyquist Velocity [m/s] or Maximum unambiguous range [m]
    wavelength : float
        Radar wavelength [m]
    """
    return (speed_of_light * wavelength / 8.) / np.asarray(varin)

    ######################
    #  Mobile Platforms  #
    ######################


def vel_shift(ground_speed, psi):
    """
    Adjusted Doppler velocity [m/s] from a mobile platform.
    Shift in Doppler velocity from mobile perspective.

    Jorgensen (1983), Eqn 2

    Parameters
    ----------
    ground_speed : float or array
        Gound speed [m/s]
    psi : float or array
        Angle between actual azimuth and fore/aft angle [deg]

    Notes
    -----
    In the case of a mobile platform (such as the NOAA P-3 aircraft, the
      Doppler velocity must be adjusted for movement of the scanning platform.

    The fore/aft angle is defined as the angle fore or aft from a plane
      normal to the direction of motion
    """
    len_gs = len(np.asarray(ground_speed))
    len_psi = len(np.asarray(psi))
#    if  len_gs != len_psi and len_gs != 1 and len_psi != 1:
    return np.asarray(ground_speed) * np.cos(np.deg2rad(psi))


def vmax_dual(wavelength, prf1, prf2):
    """Doppler velocity [m/s] from dual PRF scheme radar (+ or -).

    From Jorgensen (1983), Eqn 2

    Parameters
    ----------
    wavelength : float
        Radar wavelength [m]
    prf1 : float
        First Pulse repetition frequency [Hz]
    prf2 : float
        Second Pulse repetition frequency [Hz]

    Notes
    -----
    In the case of a mobile platform (such as the NOAA P-3 aircraft, the
    Doppler velocity must be adjusted for movement of the scanning platform.

    The fore/aft angle is defined as the angle fore or aft from a plane
    normal to the direction of motion
    """
    return (wavelength / (4 * ((1. / prf1) - (1. / prf2))))


def dual_nyquist(wavelength, prf1, prf2):
    """Wrapper for vmax_dual."""
    return vmax_dual(wavelength, prf1, prf2)


##############
#  Geometry  #
##############


def r_effective(dndh=-39e-6):
    """
    Effective radius [m] calculation.

    Rinehart (1997), Eqn 3.9, solved for R'

    Parameters
    ----------
    dndh : float
        Refraction [N x10^-6/km]

    Notes
    -----
    Effective radius of earth given a refraction.  If no refraction is given
    a "standard atmosphere" is assumed, the valued needed to have straight
    radar rays.
    """
    # Convert earth's radius to km for common dN/dH values and then
    # multiply by 1000 to return radius in meters
    return (1. / ((1/(earth_radius/1000.)) + (dndh))) * 1000.


def half_power_radius(range, bwhalf):
    """
    Half-power radius [m].

    Battan (1973),

    Parameters
    ----------
    range : float or array
        Range [m]
    bwhalf : float
        Half-power beam width [degrees]
    """
    # Convert earth's radius to km for common dN/dH values and then
    # multiply by 1000 to return radius in meters
    return (np.asarray(range) * np.deg2rad(bwhalf)) / 2.


def ray_height(range, elev, h0, reff=r43):
    """
    Center of radar beam height [m] calculation.

    Rinehart (1997), Eqn 3.12, Bech et al. (2003) Eqn 3

    Parameters
    ----------
    range : float or array
        Range from radar to point of interest [m]
    elev : float
        Elevation angle of radar beam [deg]
    h0 : float
        Height of radar antenna [m]
    reff : float
        Effective radius

    Notes
    -----
    If no Effective radius is given a "standard atmosphere" is assumed,
    the 4/3 approximation.

    Bech et al. (2003) use a factor ke that is the ratio of earth's radius
    to the effective radius (see r_effective function) and Eqn 4 in B03
    """
    # Convert earth's radius to km for common dN/dH values and then
    # multiply by 1000 to return radius in meters
    term1 = (np.sqrt(np.asarray(range)**2 + reff**2 +
             2 * np.asarray(range) * reff * np.sin(np.deg2rad(elev))))
    h = term1 - reff + h0
    return h


def sample_vol_ideal(r, bw_h, bw_v, pulse_length):
    """
    Idealized Sample volume [m^3] assuming all power in half-power beamwidths.

    From Rinehart (1997), Eqn 5.2

    Parameters
    ----------
    r : float or array
        Distance to sample volume from radar [m]
    bw_h : float
        Horizontal beamwidth [deg]
    bw_v : float
        Vertical beamwidth deg]
    pulse_length : float
        Pulse length [m]

    Notes
    -----
    This form assumes all transmitted energy is in the half-power beamwidths.
    A more realistic solution is found in the sample_vol_gauss function
    """
    return (np.pi * (np.asarray(r) * np.deg2rad(bw_h)/2.) * (np.asarray(r) *
            np.deg2rad(bw_v)/2.) * (pulse_length/2.))


def sample_vol_gauss(range, bw_h, bw_v, pulse_length):
    """
    Sample volume [m^3] assuming transmitted energy in Gaussian beam shape.

    From Rinehart (1997), Eqn 5.4

    Parameters
    ----------
    range  : float or array
        Distance to sample volume from radar [m]
    bw_h : float
        Horizontal beamwidth [deg]
    bw_v : float
        Vertical beamwidth deg]
    pulse_length : float
        Pulse length [m]

    Notes
    -----
    This form assumes a Gaussian beam shape for transmitted energy and is more
    realistic than the sample_vol_ideal.  Derived by Probert-Jones (1962).
    """
    numer = (np.pi * np.asarray(range)**2 * np.deg2rad(bw_h) *
             np.deg2rad(bw_v) * pulse_length)
    denom = 16. * np.log(2)

    SVol = numer / denom
    return SVol


def range_correct(range, height, elev):
    """
    A corrected range [m] from radar that takes into account the "loss" of
    ground distance because of the radar elevation angle.  This is a
    cumulative effect at each gate along the ray.

    From CSU Radar Meteorology AT 741 Notes

    Parameters
    ----------
    range  : float or array
        Distance to sample volume from radar [m]
    height : float
        Height of the center of radar volume [m]
    elev : float
        Elevation angle [deg]

    Notes
    -----
    This function requires that an array be passed!  If you need just one
       point create a 2 element array with a begin point.
    This is now set up to only accept a 1D array I believe.  May need to
       fix this in the future.
    """
    # Calculate the change in height along the ray
    dh1 = hheight[1:] - hheight[:-1]
    # Add the 0th place in the ray at the beginning
    dh2 = np.insert(dh1, 0, hheight[0])

    # Calculate the change in distance at each gate
    a90r = np.pi/2.  # 90 degrees in radians
    dr = dh2 / (np.tan(a90r - np.deg2rad(elev)))

    # Now calculate the corrected range at each gate
    rnew = np.asarray(range) - np.cumsum(dr)
    return rnew


def beam_block_frac(terrain, beam_height, a):
    """Partial beam blockage fraction. Unitless.

    From Bech et al. (2003), Eqn 2 and Appendix

    Parameters
    ----------
    terrain : float
        Terrain height [m]
    beam_height : float
        Beam height [m]
    a : float
        Half power beam radius [m]

    Notes
    -----
    This procedure uses a simplified interception function where no vertical
      gradient of refractivity is considered.  Other algorithms treat this
      more thoroughly.  However, this is accurate in most cases other than
      the super-refractive case.

    See the the half_power_radius function to calculate variable a

    The heights must be the same units!
    """

    # First find the difference between the terrain and height of
    # radar beam (Bech et al. (2003), Fig.3)
    y = terrain - beam_height

    numer = ((y * np.sqrt(a**2 - y**2)) + (a**2 * np.arcsin(y/a)) +
             (np.pi * a**2 / 2.))

    denom = np.pi * a**2
    return numer / denom

###############
#  Variables  #
###############


def reflectivity(power_t, gain, pulse_width, wavelength, bw_h, bw_v,
                 aloss, rloss, power_return, range, dielectric=0.93):
    """
    Radar reflectivity [mm^6/m^3].

    From Rinehart (1993), Eqn 5.17 (See Eqn 5.14-5.16 also)

    Parameters
    ----------
    power_t : float
        Transmitted power [W]
    gain : float
        Antenna Gain [dB]
    pulse_width : float
        Pulse Width [s]
    wavelength : float
        Radar wavelength [m]
    bwidth : float
        Antenna beamwidth [degrees]
    aloss : float
        Antenna/waveguide/coupler loss [dB]
    rloss : float
        Receiver loss [dB]
    dielectric : float
        Dielectric factor [unitless]
    power_return : float
        Returned power [W]
    range : float or array
        Range to target [m]

    Notes
    -----
    This routine calls the radar_constant function.

    The default is for a dielectric factor value for water.  This can be
    changed by the user, e.g. K=0.208 for particle sizes of equivalent melted
    diameters or K=0.176 for particle sizes of equivalent ice spheres.
    """
    # Call the radar constant function
    C1 = radar_constant(power_t, gain, pulse_width, wavelength, bw_h, bw_v,
                        aloss, rloss)
    return power_return * np.asarray(range)**2 / (C1 * dielectric**2)


def radial_velocity(frequency, wavelength):
    """
    Radial velocity [m/s].

    From Rinehart (1997), Eqn 6.6

    Parameters
    ----------
    frequency : float or array
        Frequency shift [Hz]
    wavelength : float or array
        Radar wavelength [m]

    Notes
    -----
    If arrays are used, either for one okay, or both must be the same length.
    """
    return np.asarray(frequency) * np.asarray(wavelength) / 2.


def cdr(refl_parallel, refl_orthogonal):
    """
    Circular depolarization ratio [dB].

    From Rinehart (1997), Eqn 10.2

    Parameters
    ----------
    refl_parallel : float or array
        Reflectivity in the parallel channel [mm^6/m^3]
    refl_orthogonal : float or array
        Reflectivity in the orthogonal channel [mm^6/m^3]

    Notes
    -----
    Ensure that both powers have the same units!

    Radars that transmit right-hand circular polarization and receive and
       receive both left- and right-hand circular polarization (using two
       antennas) and acquiring the same pulse.
    The parallel (orthogonal) component refers to the same
      (opposite) polarization as transmitted. Non-spericity of hydrometeors
      may be detected (inf long, thin scatterers have CDR = 0 dB, while perfect
     spheres have CDR = -infinity

    Can also use power measurements instead of reflectivity.

    If arrays are used, either for one okay, or both must be the same length.
    """
    return 10. * np.log10(
        np.asarray(refl_parallel)/np.asarray(refl_orthogonal))


def ldr(z_h, z_v):
    """
    Linear depolarization ratio [dB].

    From Rinehart (1997), Eqn 10.3

    Parameters
    ----------
    z_h : float or array
        Horizontal reflectivity [mm^6/m^3]
    z_v : float or array
        Vertical reflectivity [mm^6/m^3]

    Notes
    -----
    Ensure that both powers have the same units!

    Uses both polarizations in a dual-pol radar from a single pulse.

    Perfect spheres yield LDR => -infinity (though antenna limitations limit
    LDR values to -40 dB for small spheres.
    Long, thin targets, LDR => 0

    Typical values in the range -15 > LDR > -35 dB

    If arrays are used, either for one okay, or both must be the same length.
    """
    return 10. * np.log10(np.asarray(z_h) / np.asarray(z_v))


def zdr(z_h, z_v):
    """
    Differential reflectivity [dB].

    From Rinehart (1997), Eqn 10.3 and Seliga and Bringi (1976)

    Parameters
    ----------
    z_h : float or array
        Horizontal reflectivity [mm^6/m^3]
    z_v : float or array
        Vertical reflectivity [mm^6/m^3]

    Notes
    -----
    Ensure that both powers have the same units!

    Alternating horizontally and linearly polarized pulses are averaged.

    Notes
    -----
    If arrays are used, either for one okay, or both must be the same length.
    """
    return 10. * np.log10(np.asarray(z_h) / np.asarray(z_v))


def zdp(z_h, z_v):
    """
    Reflectivity difference [dB].

    From Rinehart (1997), Eqn 10.3

    Parameters
    ----------
    z_h : float
        Horizontal reflectivity [mm^6/m^3]
    z_v : float
        Horizontal reflectivity [mm^6/m^3]

    Notes
    -----
    Ensure that both powers have the same units!

    Alternating horizontally and linearly polarized pulses are averaged.
    """
    zh = np.atleast_1d(z_h)
    zv = np.atleast_1d(z_v)
    if len(zh) != len(zv):
        raise ValueError('Input variables must be same length')
        return

    zdp = np.full_like(zh, np.nan)
    good = np.where(zh > zv)
    zdp[good] = 10. * np.log10(zh[good] - zv[good])
    return zdp


def hdr(dbz_h, zdr):
    """
    Differential reflectivity [dB] hail signature.

    From Aydin et al. (1986), Eqns 4-5

    Parameters
    ----------
    dbz_h : float or array
        Horizontal reflectivity [dBZ]
    zdr : float or array
        Differential reflectivity [dBZ]

    Notes
    -----
    Ensure that both powers have the same units!

    Positive HDR and strong gradients at edges signify ice. The larger HDR,
       the greater likelihood that ice is present.

    Considerations for this equation (see paper for more details):
       1) Standar error of disdrometer data allowed for
       2) Drop oscillation accounted for based on 50% model of
          Seliga et al (1984)
       3) Lower (27) bound chose to provide constant Zh ref level
       4) Upper cutoff of 60 (may be too low)

    Picca and Ryzhkof (2012) mention that this does not take into account
    the hail melting process.  So use at your own risk!
    """
    zdr = np.atleast_1d(zdr)
    # Set the f(zdr) based upon observations
    f = np.full_like(zdr, np.nan)
    negind = np.where(zdr <= 0)
    lowind = np.where((zdr > 0) & (zdr <= 1.74))
    highind = np.where(zdr > 1.74)
    f[negind] = 27.
    f[lowind] = 19. * zdr[lowind] + 27.
    f[highind] = 60.
    # Calculate HDR
    return np.asarray(dbz_h) - f

#################
#  Attenuation  #
#################


def abs_coeff(diam, wavelength, m):
    """
    Absorption coefficient of a spherical particle. Unitless.

    Doviak and Zrnic (1993), Eqn 3.14a or Battan (1973), Eqn 6.6

    Parameters
    ----------
    diam : float or array
        Particle diameter [m]
    wavelength : float
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
    Km_im = np.imag(-1 * k_complex(m))
    Qa = (np.pi**2 * np.asarray(diam)**3 / wavelength) * Km_im
    return Qa


def scat_coeff(diam, wavelength, m):
    """
    Scattering coefficient of a spherical particle. Unitless.

    Doviak and Zrnic (1993), Eqn 3.14b or Battan (1973), Eqn 6.4

    Parameters
    ----------
    diam: float or array
        Particle diameter [m]
    wavelength: float
        Radar wavelength [m]
    m: float
        Complex refractive index [unitless]
    """
    Km_abs = np.absolute(k_complex(m))
    Qs = 2 * np.pi**5 * np.asarray(diam)**6 / (3 * wavelength**4) * (Km_abs**2)
    return Qs


def ext_coeff(diam, wavelength, m):
    """
    Extinction coefficient of a spherical particle. Unitless.

    Doviak and Zrnic (1993), Eqn 3.14b or Battan (1973), Eqn 6.5

    Parameters
    ----------
    diam : float or array
        Particle diameter [m]
    wavelength : float
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
    Qa = abs_coeff(np.asarray(diam), wavelength, m)
    Qs = scat_coeff(np.asarray(diam), wavelength, m)
    Qe = Qa + Qs
    return Qe


def k_complex(m):
    """
    Parameters
    ----------
    m : float
        Complex refractive index [unitless]
    """
    return (m**2 - 1) / (m**2 + 2)
