"""
Timothy James Lang
tjlangco@gmail.com

Last Updated 06 March 2015 (Python)
Last Updated 26 July 2005 (IDL)
"""
import numpy as np
from warnings import warn

#NEXT THREE DATA STATEMENTS CONTAIN THE FIR FILTER ORDER, GAIN, AND COEFICIENTS
#The specification of FIR filter coeficients is set for gate spacing of 150 meters.
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
    """
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
    xx = np.zeros(31) #Range window used to calculate specific differential phase
    yy = np.zeros(31) #Phase window used to calculate specific differential phase
    half = 5 #Half the window size for calculating standard deviation of phase
    half1 = 10 #Half the window size for the FIR filtering
    tmp = np.zeros(11) #Dummy variable for the standard deviation routine
    y = np.zeros(length) + bad #Dummy variable to store filtered phase
    z = 1.0 * dp #Dummy variable to store unprocessed phase
    fill = np.zeros(21) #Used to fill small gaps in phase data (phase)
    xf = np.zeros(21) #Used to fill small gaps in phase data (range)

    #Calculate standard deviation of phidp
    for i in xrange(length):
        ind = 0
        for k in np.int32(np.linspace(i-half, i+half, 2*half+1)):
            if k >= 0 and k < length:
                if dp[k] >= -180.0:
                    tmp[k-i+half] = dp[k]
                    ind += 1
        if ind > half:
            sd_lin[i] = np.std(tmp[np.arange(ind)])

    #------------- MAIN LOOP of Phidp Adaptive Filtering --------------------
    for mloop in np.arange(nfilter):
    #FIR FILTER SECTION
        for i in np.int32(np.linspace(half1, length-half1-1, length-2*half1)):
            it = 0
            if sd_lin[i] <= thsd[i] and z[i] >= -180.0:
                #Check how many points there are and prepare fill arrays
                for j in np.arange(fir3order+1):
                    if sd_lin[i-fir3order/2+j] <= thsd[i-fir3order/2+j] and\
                       z[i-fir3order/2+j] >= -180.0:
                        fill[it] = z[i-fir3order/2+j]
                        xf[it] = rng[i-fir3order/2+j]
                        it += 1
                #Enough points to fill bad areas
                if it > 16:
                    xfin = np.arange(it)
                    result = np.polyfit(xf[xfin], fill[xfin], 1)
                    acc=0.0
                    for j in np.arange(fir3order+1):
                        if sd_lin[i-fir3order/2+j] > thsd[i-fir3order/2+j] or\
                           z[i-fir3order/2+j] < -180.0:
                            z[i-fir3order/2+j] = result[0] * \
                                  rng[i-fir3order/2+j] + result[1]
                        acc += fir3coef[j] * z[i-fir3order/2+j]
                    y[i] = acc * fir3gain
                #Not enough points
                else:
                    y[i] = bad
            #Gate is bad
            else:
                y[i] = bad
    dp_lin = 1.0 * y #Store filtered PHIDP
    #*****************END LOOP for Phidp Adaptive Filtering*******************
    
    #CALCULATE KDP
    for i in np.arange(length):
        if sd_lin[i] <= thsd[i] and sd_lin[i] > 0:
            #Default value for nadp is 10, but varies based on Zh
            nadp = 10
            if i > 15 and i < length-15:
                if dz[i] < 35.0:
                    nadp = 30
                if dz[i] >= 35.0 and dz[i] < 45.0:
                    nadp = 20
                xxi=0
                for jj in np.arange(nadp+1):
                    if sd_lin[i-nadp/2+jj] <= thsd[i-nadp/2+jj] and\
                       dp_lin[i-nadp/2+jj] > -180.0:
                        xx[xxi] = rng[i-nadp/2+jj]
                        yy[xxi] = dp_lin[i-nadp/2+jj]
                        xxi += 1
                #Improved Kdp base on LSE fit to Adap filt Phidp
                if np.float(xxi)/np.float(nadp) >= 0.8:
                    xxin = np.arange(xxi)
                    result = np.polyfit(xx[xxin], yy[xxin], 1)
                    kd_lin[i] = 0.5 * result[0]
    #*******************END KDP CALCULATION******************************
    
    return kd_lin, dp_lin, sd_lin


