"""
Brenda Dolan CSU, October 2018
bdolan@atmos.colostate.edu

This is a melting layer detection algorithm. It is basically a three category
HID where it is finding the melting layer vs. anything else (other).

There are two different methodologies based on if it is a ppi or rhi scan.

Porting over Elizabeth Thompsons's HID code from IDL.
(Based on Thompson et al. 2014)

Thompson, E. J., Rutledge, S. A., Dolan, B., Chandrasekar, V., & Cheong, B. L. (2014). 
A dual-polarization radar hydrometeor classification algorithm for winter precipitation. 
Journal of Atmospheric and Oceanic Technology, 31(7), 1457-1481.

(Also based on Dolan and Rutledge, 2009)

"""

from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import scipy.stats as stats


from csu_fhc_ml1 import csu_fhc_cold_newml1


def melting_layer(dz=None,zdr=None,ldr=None,kdp=None,rho=None,sn=None,heights=None,scan_type = 'ppi',
                verbose=False,plot_flag=False,band='S',fdir='./',minRH = 0.5,expected_ML=None,nsect=36,sn_thresh=5,azimuths=None):
            
#The first step ist to run a FHC that determines 'wet snow' from 'other'
    scores = csu_fhc_cold_newml1(dz=dz, zdr=zdr, rho=rho, kdp=kdp, use_temp=False, band='S',method='linear',
                            verbose=True,fdir='/Users/bdolan/scratch/WINTER_HCA/beta_function_parameters/')

    fh = np.argmax(scores, axis=0) + 1
    
#Now this works a lot better if the SNR is included.
    if sn is not None:
        print('QCing using SN<',sn_thresh,' to remove edge WS classifications')
        fh[sn<sn_thresh] = 0
    else:
        if verbose==True:
            print('Wet snow detection will be imporved by including the SN field.')

#Now for the ML algorithm Rho is the most important variable. So if the rho was not valid, and where rho is less than the threshold. 
    rhfill = rho.filled(fill_value = np.nan)
    whbad = np.where(np.isnan(rhfill))
    fh[whbad] = -1
    print('Qcing using RH<',minRH)
    fh[rho<=minRH] = -1
        
#In this case, fh =1 is 'other', fh = 2 is 'Wet snow'
# Now use the returned Wet snow field to look at melting layer statistics and try to refine
# The location of the ML. There are 2 different methods based on if the scan is an rhi or ppi.
    if heights is None:
        print('No heights specified. Returning')
        Return            
    if scan_type == 'rhi':
        print('scan type is RHI')
        meltlev,mean_melt = get_ml_rhi(fh,dz,heights,expected_ML,verbose=verbose)
        print('RADAR MELTING IS: ',mean_melt)
    else:
        print('scan type is PPI')
        meltlev,mean_melt = get_ml_ppi(fh,dz,heights,expected_ML,azimuths=azimuths,nsect=nsect,verbose=verbose)
        print('RADAR MELTING IS: ',np.mean(mean_melt))
    #Lastly, return the HID returned from the WS / other
#     fh[meltlev==0] = 0
#     fh[meltlev==2] = 0
    
    

    return meltlev,mean_melt,fh,scores

def get_ml_rhi(fh,dz,height,expected_ML,verbose=False):
#THis will return two things. The first is an array the same size as FH data identifying the melting region, the cold and warm regions
# meltlev = 0   Warm
# meltlev = 1   Melting Layer
# meltlev = 2   Cold layer

#It will also return the height of the melting layer.

    if expected_ML is not None:
        if verbose == True:
            print('Using an expected ML of {m} to narrow the results'.format(m=expected_ML))
        wh_ht = np.ma.logical_and(height<=expected_ML+1.5,height>=expected_ML-1.5)
        wh_ws = np.ma.logical_and(wh_ht,fh ==2)
        wh_ML = np.ma.where(wh_ws)
    else:
        wh_ML = np.ma.where(fh==2)
# find total number of ML pixels within these "proper" limits of height and range
    ML_z = height[wh_ML]
# define 80th and 20th percentile marks for top and bottom of ML ~ giangrande et al. 2008 / boodoo et al. 2010 methodology
    ML_num = len(ML_z)

###### you are going to have some error messages if it doesn't find a melting layer - it needs an array to do these statistics
# this algorithm is designed for stratiform precipitation, so if it doesn't have a ML you really should be using a different algorithm
# because there are no convective precip categories like graupel, hail, etc.
    if len(ML_z) > 10:
        distfreq,eg = np.histogram(ML_z)
        maxfreq=np.max(distfreq)
        ML_mode,ML_count = stats.mode(ML_z)
        ML_sdevZ = np.std(ML_z)
        ML_medianZ = np.median(ML_z)
        ML_varZ = np.var(ML_z)
        ML_skewZ = stats.skew(ML_z)
        ML_kurtZ = stats.kurtosis(ML_z)
        try:
            ML_80Z = np.percentile(ML_z,80)
        except:
            if verbose == True:
                print('Not enough points for 80 percent')
            else:
                pass
        try:
            ML_20Z = np.percentile(ML_z,20)
        except:
            if verbose == True:
                print('Not enough points for 20 percent')
            else:
                pass
        if verbose == True:
            print("ML Stats:")
            print("mode: {m:.2f}".format(m=ML_mode[0]))
            print("Median:{o:.2f}".format(o=ML_medianZ))
            print("SDEV: {s:.2f}".format(s=ML_sdevZ))
            print("Var: {v:.2f}".format(v=ML_varZ))
            print("Kurt: {k:.2f}".format(k=ML_kurtZ))
            print("Skew: {w:.2f}".format(w=ML_skewZ))
            print("80per: {e:.2f}".format(e=ML_80Z))
            print("20per: {t:.2f}".format(t=ML_20Z))
        else:
            print("Melting level height:{m:.2f}".format(m=ML_medianZ))
    else:
        print("not enough points for ML detection. Using expected: {m}".format(m=expected_ML))
        ML_medianZ = expected_ML
    
    meltlev = np.zeros_like(fh)+1
    
    wh_warm = np.where(height  < ML_medianZ)
    wh_cold = np.where(height > ML_medianZ)

    meltlev[wh_warm] = 0
    meltlev[wh_cold] = 2
    
    if verbose == True:
        print("ML Stats:")
        print("mode: {m:.2f}".format(m=ML_mode[0]))
        print("Median:{o:.2f}".format(o=ML_medianZ))
        print("SDEV: {s:.2f}".format(s=ML_sdevZ))
        print("Var: {v:.2f}".format(v=ML_varZ))
        print("Kurt: {k:.2f}".format(k=ML_kurtZ))
        print("Skew: {w:.2f}".format(w=ML_skewZ))
        print("80per: {e:.2f}".format(e=ML_80Z))
        print("20per: {t:.2f}".format(t=ML_20Z))
    else:
        print("Melting level height:{m:.2f}".format(m=ML_medianZ))
    
    #lastly, return the hid associated with WS and Other
    
    return meltlev,ML_medianZ
    
def get_ml_ppi(fh,dz,height,expected_ML,nsect=36,verbose=False,azimuths=None):

# This is defining variable ML height for each azimuth sector!
# 36 sectors with 10 deg intervals
#THis will return two things. The first is an array the same size as FH data identifying the melting region, the cold and warm regions
# meltlev = 0   Warm
# meltlev = 1   Melting Layer
# meltlev = 2   Cold layer

#It will also return the height of the melting layer.
    if azimuths is not None:
        az3d = np.repeat(azimuths[...,np.newaxis],np.shape(height[0,:]),axis=1)
        nsec=36
        if verbose == True:
            print('Number of sections:{n}'.format(n=nsec))
        nint = 36
        az_sector = np.linspace(0,360,nint)
        ML_mode = np.zeros([nint-1])
        ML_sdevZ = np.zeros([nint-1])
        ML_medianZ = np.zeros([nint-1])
        ML_varZ = np.zeros([nint-1])
        ML_skewZ = np.zeros([nint-1])
        ML_kurtZ = np.zeros([nint-1])
        ML_80Z = np.zeros([nint-1])
        ML_20Z = np.zeros([nint-1])
        meltlev_bottom_az = np.zeros_like(fh).astype(float)
        meltlev_top_az = np.zeros_like(fh).astype(float)
        meltlev_az = np.zeros_like(fh).astype(float)
        meltlev= np.zeros_like(fh).astype(float)+1.


        for i,a in enumerate(az_sector[:-1]):
            wh_az = np.logical_and(az3d>=a,az3d<az_sector[i+1])
            whsec = np.squeeze(np.ma.where(wh_az))
            wh_sector1 = np.ma.where(np.logical_and(wh_az,fh ==2))
            wh_ws = np.ma.where(fh[wh_sector1]==2)

            ML_sector = height[wh_sector1]
            npts = len(np.ravel(whsec))
            #print(npts)
            ML_num_sector = len(ML_sector)
            per_ML = np.float(ML_num_sector)/npts*100.
            #print('ML number in sector',ML_num_sector,'percent of sector',per_ML)
            #print('Mean dz in ML',np.nanmean(dz[wh_sector1]))
            if (per_ML > 0.05) and (np.nanmean(dz[wh_sector1])> 15):
                #print 'made it to if'
                #print i,a
                MLcheck = np.median(ML_sector[wh_ws])
                print (MLcheck,expected_ML)
                if np.abs(MLcheck-expected_ML)>2.:
                    print('Melting layer is way off. Double check sounding and radar data or consider different HCA.')
                else:
                    ML_mode[i],ML_count = stats.mode(ML_sector[wh_ws])
                    ML_sdevZ[i] = np.std(ML_sector[wh_ws])
                    ML_medianZ[i] = np.median(ML_sector[wh_ws])
                    ML_varZ[i] = np.var(ML_sector[wh_ws])
                    ML_skewZ[i] = stats.skew(ML_sector[wh_ws])
                    ML_kurtZ[i] = stats.kurtosis(ML_sector[wh_ws])
                    ML_80Z[i] = np.percentile(ML_sector[wh_ws],90)
                    ML_20Z[i] = np.percentile(ML_sector[wh_ws],10)
            else:
                if verbose == True:
                    print("Too few points for ML. Consider warm HCA or completely winter case.")
            whoutside = np.squeeze(np.where(np.logical_or(ML_sector<=(ML_20Z[i]),ML_sector>=(ML_80Z[i]))))
            meltlev_bottom_az[whsec]=ML_20Z[i]
            meltlev_top_az[whsec]=ML_80Z[i]
            meltlev_az[whsec]=ML_medianZ[i]
            if verbose == True:
                print('Melting layer for {a}-{b}: {m}'.format(a=a,b=az_sector[i+1],m=ML_medianZ[i]))

        wh_warm = np.where(height  <= meltlev_az)
        wh_cold = np.where(height  >= meltlev_az)
        meltlev[wh_cold] = 2
        meltlev[wh_warm] = 0
#        meltlev[wh_sector1]=meltlevlim
    else:
        print('No azimuths specified. Using RHI version.')
        meltlev,ML_medianZ = get_ml_rhi(fh,height,expected_ML,verbose=False)
    return meltlev,ML_medianZ