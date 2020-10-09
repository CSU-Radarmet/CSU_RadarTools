"""
csu_fhc_winter.py

Brenda Dolan CSU, October 2018
bdolan@atmos.colostate.edu

Porting over Elizabeth Thompsons's HID code from IDL
(Based on Thompson et al. 2014)
(Also based on Dolan and Rutledge, 2009)

Thompson, E. J., Rutledge, S. A., Dolan, B., Chandrasekar, V.,
& Cheong, B. L. (2014). A dual-polarization radar hydrometeor classification
algorithm for winter precipitation. Journal of Atmospheric and Oceanic
Technology, 31(7), 1457-1481.

Cateories:
0  = Unclassified
1  = Ice Crystals
2  = Plates
3  = Dendrites
4  = Aggregates
5  = Wet Snow
6  = Frozen precip
7  = Rain

"""

from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from .beta_functions import get_mbf_sets_winter, CSV_DIR
from .calc_kdp_ray_fir import hid_beta_f
from .csu_fhc_melt import melting_layer
from .csu_fhc_winter import csu_fhc_winter


def run_winter(dz=None, zdr=None, rho=None, kdp=None, ldr=None, sn=None,
               T=None, use_temp=False, band='S', method='linear', sn_thresh=5,
               expected_ML=4.0, nsect=36, return_scores=False, azimuths=None,
               verbose=True, minRH=0.5, scan_type='ppi', heights=None,
               fdir=CSV_DIR):

    # The first step is to find the wet snow and the melting layer.

    meltlev, melt_z, fh, scores_ML = melting_layer(
        dz=dz, zdr=zdr, kdp=kdp, rho=rho, sn=sn, heights=heights,
        scan_type=scan_type, verbose=verbose, band=band,
        fdir=fdir, azimuths=azimuths, expected_ML=expected_ML, minRH=minRH,
        sn_thresh=sn_thresh, nsect=36)

    # Step 2 is to run the warm layer HID.
    # Using csu_fhc_winter and warm == True

    scores_warm = csu_fhc_winter(
        dz=dz, zdr=zdr, rho=None, kdp=kdp, use_temp=True, T=T, band='C',
        warm=True, verbose=True, fdir=fdir)
    fhwarm = np.argmax(scores_warm, axis=0) + 1
    # Now reset the fhwarm values to correspond to
    # 6 - Frozen
    # 7 - Rain
    fhwarm[fhwarm == 1] = 6
    fhwarm[fhwarm == 2] = 7

    scores_cold = csu_fhc_winter(dz=dz, zdr=zdr, rho=None, kdp=kdp,
                                 use_temp=True, T=T, band='C', warm=False,
                                 verbose=True, fdir=fdir)

    fhcold = np.argmax(scores_cold, axis=0) + 1

    # Now reset the fh values to correspond to
    # 5 - Wet Snow
    whgd = np.where(fh == 2)
    fh[whgd] = 5

    winter_hca = np.zeros_like(fh) - 1

    # Now combine them using the melting level idea.
    winter_hca[meltlev == 0] = fhwarm[meltlev == 0]
    winter_hca[meltlev == 2] = fhcold[meltlev == 2]
    # whbad = np.where(fh == -1)
    whmelt = np.where(fh == 5)
    winter_hca[whmelt] = 5
    winter_hca[fh == -1] = -1
    winter_hca[fh == 0] = -1

    # Filter out where there is no pol data.
    # rhfill = rho.filled(fill_value = np.nan)
    # whbad = np.where(np.isnan(rhfill))
    # winter_hca[whbad] = -1

    if return_scores:
        scores = {'ML': scores_ML, 'warm': scores_warm, 'cold': scores_cold}
        hca = {'ML': fh, 'warm': fhwarm, 'cold': fhcold}
        return winter_hca, scores, hca
    else:
        return winter_hca
