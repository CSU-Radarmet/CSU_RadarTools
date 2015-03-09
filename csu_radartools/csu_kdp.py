"""
Timothy James Lang
tjlangco@gmail.com

Last Updated 06 March 2015 (Python)
Last Updated 26 July 2005 (IDL)

Python v1.0

"""
import numpy as np
from numpy import linalg
from warnings import warn
#import time

VERSION = '1.0'
#NEXT THREE DATA STATEMENTS CONTAIN THE FIR FILTER ORDER, GAIN, AND COEFICIENTS
#The specification of FIR filter coeficients is set for gate spacing of
#150 meters.
fir3order = 20
fir3gain = 1.044222
fir3coef = np.array([1.625807356e-2, 2.230852545e-2, 2.896372364e-2,
                     3.595993808e-2, 4.298744446e-2, 4.971005447e-2,
                     5.578764970e-2, 6.089991897e-2, 6.476934523e-2,
                     6.718151185e-2, 6.800100000e-2, 6.718151185e-2,
                     6.476934523e-2, 6.089991897e-2, 5.578764970e-2,
                     4.971005447e-2, 4.298744446e-2, 3.595993808e-2,
                     2.896372364e-2, 2.230852545e-2, 1.625807356e-2])

def calc_kdp_bringi(dp=None, dz=None, rng=None, thsd=12, nfilter=1, bad=-32768):
    """
    Overview
    --------
    This is a decade+ old algorithm that uses a 21-pt FIR filter to process
    differential phase and extract specific differential phase.
    It works on individual rays of polarimetric radar data. It is based on code
    provided by V. N. Bringi and Yanting Wang of CSU Electrical Engineering.
    It assumes differential phase has been unfolded already.
    
    Steps
    -----
    1. Standard deviation of differential phase is calculated and used to QC the 
       phase data.
    2. Differential phase is filtered. 
    3. Specific differential phase is calculated by consulting reflectivity. 
       As reflectivity declines progressively more and more gates are needed
       in the window used to fit a line to the filtered phase. Specific
       differential phase is half the slope of that line.

    Reference
    ---------
    Timothy J. Lang, David A. Ahijevych, Stephen W. Nesbitt, Richard E. Carbone, 
    Steven A. Rutledge, and Robert Cifelli, 2007: Radar-Observed Characteristics of 
    Precipitating Systems during NAME 2004. J. Climate, 20, 1713–1733. 
    doi: http://dx.doi.org/10.1175/JCLI4082.1
    
    Arguments
    ---------
    dp = Differential phase (deg, vector)
    dz = Reflectivity (dBZ, vector)
    rng = Range (km, vector)
    thsd = Threshold for standard deviation of differential phase, above which 
           the data are not considered when filtering or calculating specific
           differential phase. The user can specify a vector of spatially
           varying thresholds instead (e.g., vary by range).
    nfilter = Number of times to apply the FIR filter
    bad = Value for bad/missing data
    
    Returns
    -------
    kd_lin = Specific differential phase (deg/km, vector)
    dp_lin = Filtered differential phase (deg, vector)
    sd_lin = Standard deviation of differential phase (deg, vector)
    
    To Do
    -----
    1. Performance improvements
    2. Make object-oriented
    """
    #begin_time = time.time()
    #Quick check for all the variables. Used keywords so order doesn't matter.
    if dp is None or dz is None or rng is None:
        warn('Missing needed variables (dp, dz, rh, and/or rng), failing ...')
        return
    #Define needed variables
    kd_lin = np.zeros_like(rng) + bad
    sd_lin = np.zeros_like(rng) + 100.0
    #User can provide a spatially varying stddev(dp) threshold
    if not hasattr(thsd, '__len__'):
        thsd = np.zeros_like(rng) + thsd
    length = len(rng)
    lin = np.arange(length)
    half = 5 #Half the window size for calculating standard deviation of phase
    half1 = fir3order / 2 #Half the window size for the FIR filtering
    y = np.zeros(length) + bad #Dummy variable to store filtered phase
    z = 1.0 * dp #Dummy variable to store un/pre-processed phase
    #print time.time() - begin_time, 'seconds since start (DEF)'

    #Calculate standard deviation of phidp
    mask = dp >= -180
    for i in lin[mask]:
        index1 = i - half
        index2 = i + half
        if index1 >= 0 and index2 < length - 1:
            yy = dp[index1:index2]
            tmp_mask = mask[index1:index2]
            if len(yy[tmp_mask]) > half:
                #Following is faster than np.std
                a = yy[tmp_mask]
                m = a.mean()
                c = a - m
                sd_lin[i] = (np.dot(c,c) / a.size)**0.5
                #sd_lin[i] = np.std(tmp[sdp_mask])
    #print time.time() - begin_time, 'seconds since start (SDP)'

    #------------- MAIN LOOP of Phidp Adaptive Filtering --------------------
    #FIR FILTER SECTION
    for mloop in np.arange(nfilter):
        mask = np.logical_and(sd_lin <= thsd, z >= -180)
        for i in lin[mask]:
            index1 = i - half1
            index2 = i + half1
            if index1 >= 0 and index2 < length - 1:
                yy = z[index1:index2+1]
                xx = rng[index1:index2+1]
                tmp_mask = mask[index1:index2+1]
                siz = len(yy[tmp_mask])
                if siz > 16:
                    if siz < 21:
                        #Following is faster than np.polyfit
                        A = np.array([xx[tmp_mask], np.ones(siz)])
                        result = linalg.lstsq(A.T, yy[tmp_mask])[0]
                        #result = np.polyfit(xx[tmp_mask], yy[tmp_mask], 1)
                        yy[~tmp_mask] = result[0] * xx[~tmp_mask] + result[1]
                    y[i] = fir3gain * np.dot(fir3coef, yy)
        z = 1.0 * y #Enables re-filtering of processed phase
    dp_lin = 1.0 * y
    #print time.time() - begin_time, 'seconds since start (FDP)'
    #*****************END LOOP for Phidp Adaptive Filtering*******************
    
    #CALCULATE KDP
    #Default value for nadp is 10, but varies based on Zh
    nadp = np.int16(0 * dz + 10)
    tmp_mask = dz < 35
    nadp[tmp_mask] = 30
    tmp_mask = np.logical_and(dz >= 35, dz < 45)
    nadp[tmp_mask] = 20
    mask = dp_lin != bad
    for i in lin[mask]:
        half_nadp = nadp[i] / 2
        index1 = i - half_nadp
        index2 = i + half_nadp + 1
        if index1 >= 0 and index2 <= length:
            tmp_mask = mask[index1:index2]
            xx = rng[index1:index2]
            siz = len(xx[tmp_mask])
            #Improved Kdp based on LSE fit to Adap filt Phidp
            if siz >= 0.8 * nadp[i]:
                yy = dp_lin[index1:index2]
                #Following is faster than np.polyfit
                A = np.array([xx[tmp_mask], np.ones(siz)])
                result = linalg.lstsq(A.T, yy[tmp_mask])[0]
                #result = np.polyfit(xx[tmp_mask], yy[tmp_mask], 1)
                kd_lin[i] = 0.5 * result[0]
    #*******************END KDP CALCULATION******************************

    #print time.time() - begin_time, 'seconds since start (KDP/Done)'
    return kd_lin, dp_lin, sd_lin


