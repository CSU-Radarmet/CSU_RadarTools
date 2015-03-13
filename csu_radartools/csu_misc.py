"""
#############
csu_misc sub-module of csu_radartools

Contacts
--------
Brenda Dolan (dolan@atmos.colostate.edu)
Paul Hein (hein@atmos.colostate.edu)
Timothy Lang (tjlangco@gmail.com)

References
----------
Lang et al. (2007; JCLIM) - Insect filter
#############
"""

import numpy as np

#Biological scatterer thresholds originally created using NAME data
DEFAULT_DZ_RANGE = [ [-100,10], [10,15], [15,20], [20,25], [25,30], [30,35] ]
DEFAULT_DR_THRESH = [1, 1.3, 1.7, 2.1, 2.5, 2.8]

def insect_filter(dz, zdr, mask=None, dz_range=DEFAULT_DZ_RANGE,
                  dr_thresh=DEFAULT_DR_THRESH):
    """
    Returns a mask that identifies potentially suspect gates due to presence of
    biological scatterers.
    """
    #Define generic mask if none provided
    if mask is None:
        mask = np.logical_or(dz <= -100, zdr < -100)
    cond = []
    for i, thresh in enumerate(dr_thresh):
        dz_subrange = dz_range[i]
        cond_dz = np.logical_and(dz > dz_subrange[0], dz <= dz_subrange[1])
        sub_cond = np.logical_and(cond_dz, zdr >= thresh)
        cond.append(sub_cond)
    store_cond = mask
    for sub_cond in cond:
        store_cond = np.logical_or(store_cond, sub_cond)
    return store_cond

def second_trip_filter_magnetron(vr, bad=-32768):
    """
    With magnetron transmitters, there is no coherent Doppler signature in 
    second trip. Thus, if velocity is missing but there is echo, that echo is 
    likely second trip.
    """
    return vr == bad

def differential_phase_filter(sdp, thresh_sdp=12):
    """
    Filter out gates with large standard deviations of differential phase.
    Note: thresh_sdp can be an array the same size as sdp.
    """
    return sdp > thresh_sdp




