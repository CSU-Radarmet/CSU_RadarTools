"""
csu_blended_rain.py

# Brody Fuchs, CSU, Oct 2014
# brfuchs@atmos.colostate.edu

# python remake of Brenda's code to calculate water and ice mass

Amendments by
Timothy Lang (tjlangco@gmail.com)
2/20/2015

# See original comments below
#***************************************************************************
#***************************************************************************
# Based on Code from A. Rowe and T. Lang from L. Nelson and L. Carey.
# Brenda Dolan
# February 27, 2012
# bdolan@atmos.colostate.edu
#***************************************************************************
#***************************************************************************

#############################
ORIGINAL FUNCTION NOTES BELOW
KEPT FOR POSTERITY ONLY
MANY POINTS OUT OF DATE
#############################
# Brody Fuchs, Oct 2014, converting the IDL code to python here
# This program essentially dublicates the fortran program
# hybridrain_mod5.f.  It determine best polarimetric 
# rainfall estimator, based on thresholds of ZH, ZDR, and KDP,
# as well as Ice fraction.  
# The input files are 2-D cdf (gridded to a specific ht).
#
# The program also does hydro ID, using the appropriate subroutines.
# Because the input cdfs are 2D (1 ht), temp info at the height of
# interest should be input for the event of interest

# fi=ice fraction
# dz=radar reflectivity
# dr=differential reflectivity
# kd=specific differntial phase
# rkd=R(kdp)
# rkddr=R(kdp,zdr)
# rzhdr=R(zh,Zdr)

# create zk "blended" rainfall (zk) array and an array to identify which
# method is used
# 1 = R(kdp,zdr)
# 2 = R(kdp)
# 3 = R(zh,zdr)
# 4 = R-Z NXRAD (total reflectivity)
# 5 = R-Z NEXRAD (rain only reflectivity)

#cdf_array=intarr(3)

### IF CALLING FROM BRORADAR, ALREADY HAVE DATA, DONT NEED TO READ IN    
# open the cdf files
#    dir = '/data3/brfuchs'
#cdffile = findfile(dir+'*.cdf',count=numfile)
#    cdffile = '%s/CHL20120606_220558_cedunf.cdf'

#    vars = ['DZ', 'DR', 'KD', 'LH', 'RH']

#    radar_data = {}
#    nc = Dataset(filename)

# READ IN NETCDF FILE AND GET DATA
#    for var in vars:
#    	radar_data[var] = nc.variables[var][:]
# create x,y arrays that will be used in the output cdf files
#    lon=xfirst+xdelta*findgen(cdf_array[0])
#    lat=yfirst+ydelta*findgen(cdf_array[1])

# BF have x and y arrays already

# Run Kyle W.'s fuzzy logic HID program

# first call the membership routine that does the fuzzification...
# use the 2004 updated version of the MBFs with the "nomix" option
#   cdf_beta_functions_kwiens,mbf_sets,/use_temp,/NoMix

# now run the HID program
# first have to define a couple of structures
# 1st structure is the radar parameter tags and array 
# that the tag corresponds to - this is described in the 
# HID subroutine cdf_fhc.pro
#
# for chill rain maps, input cdf files have only DZ,DR, and KD
#   radar_data = {DZ:dz,DR:dr,KD:kd}

# define a structure of which input variables to actually use in the HID
# don't really need this since we want to use all the variables in the HID
# process
#   fhc_vars = {Zh:1, Zdr:1, Kdp:1, LDR:0, rho_hv:0,T:1}

# CHILL rainmaps are 2D at 1 km height (CAPPI or ~ 1km height for
# SPRINT  gridppi option).  So,input a T at 1km AGL for event of
# interest (based on sounding file)
#   T=dz
#   T[*,*]=17.8

# now run the HID program over the grid and stuff the results into a
# new array (hid)
#   hid= cdf_fhc(radar_data,mbf_sets,fhc_vars,T,/use_temp,/ppi)

# here's the number of hids - useful for defining more arrays...
#   nhid=n_elements(mbf_sets.zh_set.a)
    nhid = 10

# set and initialize arrays for rr estimation

#############################
END ORIGINAL FUNCTION NOTES
#############################

#############################
ORIGINAL FUNCTION NOTES BELOW
KEPT FOR POSTERITY ONLY
MANY POINTS OUT OF DATE
#############################
###### HID rain rate calculations ########

#This is a blended polarimetric rainfall algorithm that selects 
#rainfall estimators based on a hydrometeor identification scheme.
#
#Modified by Brenda Dolan
#bdolan@atmos.colostate.edu
#June 2010
#Colorado State University

#Ported over to python by Brody Fuchs, Oct 2014

# first step is to fit the consensus of hid from all 6 inputs.  There
# are 3 categories to consider: rain, mix, ice. the
# category with the biggest number wins. 
# 1=rain(c_hid_consensus=1)
# 2=mix(c_hid_consensus=2)
# 3=ice(c_hid_consensus=3)
# 4=rain-mix tie (average of rain rr and mix rr is used)
# (c_hid_consensus=4)
# 5 rain-ice tie (rain value is used in rr calculation)
# (c_hid_consensus=5)
# 6 mix-ice tie (mix value is used in rr calculation)
# (c_hid_consensus=6)

# this is a variable that keeps track of which hid is actually
# used for the site in question, based on the relative totals 
# of the 6 values.  
#   c_hid_consensus=dz
#   c_hid_consensus[*]=0.

# here is variable to store the actual rr estimate method for hidzk
# the values go from 1-5 just like for zk; however, there is also a method
# 6 which represents the result when mix and rain HID end in a tie: 
# in this case, an average between the value calculated using the mix (r(kd))
# and rain methods (could be r(kd), r(kd,dr), r(zh) or r(zh,dr))
#############################
END ORIGINAL FUNCTION NOTES
#############################

"""

import numpy as np
import warnings
from csu_liquid_ice_mass import calc_zh_zv, calc_zdp, linearize, get_linear_fits

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

def calc_blended_rain(dz=None, zdr=None, kdp=None, ice_flag=False,
                      thresh_dz=38.0, thresh_zdr=0.5, thresh_kdp=0.3,
                      thresh_frac_ice=0.1, thresh_nexrad=53.0, a=300.0, b=1.4,
                      fit_a=None, fit_b=None, method='cr1995'):
 
    """
    This algorithm ingests polarimetric radar data and computes rain rate,
    using difference reflectivity as a guide for determining the presence of ice.
    
    Inputs:
    dz = Reflectivity
    zdr = Differential Reflectivity
    kdp = Specific Differential Phase
    fhc = Hydrometeor ID
    thresh_dz = Threshold for dz to use certain rain algorithms
    thresh_zdr = Threshold for zdr to use certain rain algorithms
    thresh_kdp = Threshold for kdp to use certain rain algorithms
    thresh_nexrad = Threshold for dz to cap Z-R    
    a, b => Z = a * R**b
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
    #Initialize, check for all vars, check for scalars
    if dz is None or kdp is None or zdr is None:
        warnings.warn('No dz, zdr, or kdp provided, failing ...')
        return
    len_flag = hasattr(dz, '__len__')
    if not len_flag:
        dz = np.array([dz])
        kdp = np.array([kdp])
        zdr = np.array([zdr])
    r_blended = np.zeros_like(dz)
    meth = np.int16(np.zeros_like(dz))

    #NEXRAD estimate
    r_dz_nexrad = calc_rain_nexrad(dz, thresh_nexrad, a=a, b=b)

    #Polarimetric estimates
    r_dz_zdr = calc_rain_z_zdr(dz, zdr)
    r_kdp = calc_rain_kdp(kdp)
    r_kdp_zdr = calc_rain_kdp_zdr(kdp, zdr)

    #calculate rain line, zdp, and ice fraction
    zhor, zvert = calc_zh_zv(dz, zdr)
    zdp = calc_zdp(zhor, zvert)
    #calculate contribution to Zh from pure rain
    if fit_a == None:
        fit_a, fit_b = get_linear_fits(method=method)
    zrain = 10.0**((fit_a * zdp + fit_b)/10.0)
    dzrain = 10.0 * np.log10(zrain)
    fi = 1.0 - (zrain / zhor)
    r_dz_rainonly = calc_rain_nexrad(dzrain, thresh_nexrad, a=a, b=b)

    #Conditions
    cond_ice = fi <= thresh_frac_ice
    cond_kdp = kdp >= thresh_kdp
    cond_dz = dz > thresh_dz
    cond_zdr = zdr >= thresh_zdr
    cond_dz_kdp = np.logical_and(cond_dz, cond_kdp)
    #Set of method choices
    cond_meth_1or2 = np.logical_and(cond_ice, cond_dz_kdp)
    cond_meth_1 = np.logical_and(cond_meth_1or2, cond_zdr)
    cond_meth_2a = np.logical_and(cond_meth_1or2, ~cond_zdr)
    cond_meth_2b = np.logical_and(~cond_ice, cond_dz_kdp)
    cond_meth_2 = np.logical_or(cond_meth_2a, cond_meth_2b)
    cond_meth_3or4 = np.logical_and(cond_ice, ~cond_dz_kdp)
    cond_meth_3 = np.logical_and(cond_meth_3or4, cond_zdr)
    cond_meth_4 = np.logical_and(cond_meth_3or4, ~cond_zdr)
    cond_meth_5 = np.logical_and(~cond_ice, ~cond_dz_kdp)
    #Assign methods
    meth[cond_meth_1] = 1
    meth[cond_meth_2] = 2
    meth[cond_meth_3] = 3
    meth[cond_meth_4] = 4
    meth[cond_meth_5] = 5
    #Assign rain rates based on methods
    r_blended[cond_meth_1] = r_kdp_zdr[cond_meth_1]
    r_blended[cond_meth_2] = r_kdp[cond_meth_2]
    r_blended[cond_meth_3] = r_dz_zdr[cond_meth_3]
    r_blended[cond_meth_4] = r_dz_nexrad[cond_meth_4]
    r_blended[cond_meth_5] = r_dz_rainonly[cond_meth_5]

    #Return based on what the user provided and what they wanted
    if not ice_flag:
        if not len_flag:
            return r_blended[0], meth[0]
        else:
            return r_blended, meth
    else:
        if not len_flag:
            return r_blended[0], meth[0], zdp[0], fi[0]
        else:
            return r_blended, meth, zdp, fi

def csu_hidro_rain(dz=None, zdr=None, kdp=None, fhc=None,
                   thresh_dz=38.0, thresh_zdr=0.5, thresh_kdp=0.3,
                   thresh_nexrad=53.0, a=300.0, b=1.4):
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
    a, b => Z = a * R**b
    
    Returns: rain rate, method
    method = 1: R(Kdp, Zdr)
    method = 2: R(Kdp)
    method = 3: R(Z, Zdr)
    method = 4: R(Z)
    
    See Bringi and Chandrasekar textbook for more information
    """
    #Initialize, check for all necessary vars, and allow scalars
    if dz is None or kdp is None or zdr is None or fhc is None:
        warnings.warn('No dz, zdr, kdp, or fhc provided, failing ...')
        return
    dz, zdr, kdp, len_flag = _check_for_array(dz, zdr, kdp)
    crr_hidzk = np.zeros_like(dz)
    crr_meth = np.int16(np.zeros_like(dz))

    #Rainfall calculations
    crr_z88t = calc_rain_nexrad(dz, a=a, b=b)
    crr_zhdr = calc_rain_z_zdr(dz, zdr)
    crr_kd = calc_rain_kdp(kdp)
    crr_kddr = calc_rain_kdp_zdr(kdp, zdr)

    #Conditions
    #fhc = 1, 2, 10 => drizzle, rain, big drops
    cond_rain_a = np.logical_or(fhc==1, fhc==2)
    cond_rain = np.logical_or(fhc==10, cond_rain_a)
    cond_kdp = kdp >= thresh_kdp
    cond_dz = dz > thresh_dz
    cond_zdr = zdr >= thresh_zdr
    cond_dz_kdp = np.logical_and(cond_dz, cond_kdp)
    cond_dz_kdp_zdr = np.logical_and(cond_dz_kdp, cond_zdr)
    cond_dz_kdp_not_zdr = np.logical_and(cond_dz_kdp, ~cond_zdr)
    cond_not_kdp_zdr = np.logical_and(~cond_kdp, cond_zdr)
    cond_not_kdp_not_zdr = np.logical_and(~cond_kdp, ~cond_zdr)
    #fhc = 5, 7, 8, 9 => wet snow, graupel (LD & HD), hail
    cond_ice_a = np.logical_or(fhc==5, fhc==7)
    cond_ice_b = np.logical_or(fhc==8, fhc==9)
    cond_ice = np.logical_or(cond_ice_a, cond_ice_b)
    cond_ice_kdp = np.logical_and(cond_ice, cond_kdp)
    #Set methods
    cond_meth_1 = np.logical_and(cond_rain, cond_dz_kdp_zdr)
    cond_meth_2 = np.logical_and(cond_rain, cond_dz_kdp_not_zdr)
    cond_meth_3 = np.logical_and(cond_rain, cond_not_kdp_zdr)
    cond_meth_4 = np.logical_and(cond_rain, cond_not_kdp_not_zdr)
    #Assign rain rates based on methods
    crr_hidzk[cond_meth_1] = crr_kddr[cond_meth_1]
    crr_hidzk[cond_meth_2] = crr_kd[cond_meth_2]
    crr_hidzk[cond_meth_3] = crr_zhdr[cond_meth_3]
    crr_hidzk[cond_meth_4] = crr_z88t[cond_meth_4]
    crr_hidzk[cond_ice_kdp] = crr_kd[cond_ice_kdp]
    #Save methods
    crr_meth[cond_meth_1] = 1
    crr_meth[cond_meth_2] = 2
    crr_meth[cond_meth_3] = 3
    crr_meth[cond_meth_4] = 4
    crr_meth[cond_ice_kdp] = 2

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
