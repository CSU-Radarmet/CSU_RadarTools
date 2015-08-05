"""
NAME:
beta_functions.py

Smaller b values means falloff is slower

PURPOSE:
Contruct Membership Beta Functions (MBFs) parameters for each
hydrometeor type and fuzzy set based on scattering simulations (Dolan and
Rutledge, 2009)

# Basic form of a Membership Beta Function is:
#
#                         1
#               ---------------------
# beta(x) =         |              |^b
#                1+ |  ((x-m)/a)^2 |
#                   |              |
#
#  Where x = input data (Zh, Zdr, etc)
#        m = center
#        a = width
#        b = slope

Anticipate having 6 input data fields and 10 fuzzy sets (HID types)

Input data fields:           HID types:           Fuzzy set:
-------------------           ----------           ----------
     Z_h                      Drizzle                  1
     Z_dr                     Rain                     2
     K_dp                     Ice Crystals             3
     LDR                      Aggregates               4
     rho_hv                   Wet Snow                 5
     Temperature              Vertical Ice             6
                              Low Density Graupel      7
                              High Density Graupel     8
                              Hail                     9
                              Big Drops               10

ORIGNAL AUTHOR:
Kyle C. Wiens

PRINCIPAL CONTACT:
Brenda Dolan
bdolan@atmos.colostate.edu

ORIGINAL DATE:
26 April 2002.

MODIFICATIONS:
28 April 2002:  Changed the fuzzy sets to include vertical ice and
get rid of low/high density snow.

15 May 2002:  Changed KDP MBF for vertical ice to be unity from -0.2
to -0.6.  Used to be from -0.6 to 0.0.  This didn't make much difference.

19 February 2003:  Changed MBFs for wet/dry graupel and Wet Snow to
conform to Liu and Chandrasekar (2000).  Changed MBF for vertical ice to
conform to Carey and Rutledge (1998).

29 January 2009:  Changed the MBFS to the theory-based S-band values. NOTE:
                  LDR VALUES WERE NOT MODIFIED.

08 February 2012: Changed the categories and the MBF values based on scattering
                  simulations. NOTE: LDR values from simulations were added. BD
10 June 2012: Adjusted MBFs for performance. BD
25 January 2015: Pythonized (TJL)
05 August 2015: Python 3 compatible (TJL)
"""

from __future__ import division
from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Beta function parameters stored in individual CSVs in separate directory
CSV_DIR = os.sep.join([os.path.dirname(__file__),
                      'beta_function_parameters'])+'/'

################################
# Helper Functions Below
################################


def _get_beta_example(var, mbf_sets, field, i):
    """Internal beta function calculator"""
    return 1.0/(1.0 + (((var - mbf_sets[field]['m'][i]) /
                mbf_sets[field]['a'][i])**2)**mbf_sets[field]['b'][i])


def plot_beta_functions(mbf_sets, band, temp_factor, pdir='./', n_types=10):
    """
    Optional plotting routine that produces basic beta function curves
    for every hydrometeor species
    mbf_sets = Membership beta function sets
    band = 'X', 'C' or 'S'
    temp_factor = Factor to modify depth of T effects
    fdir = Path to plot locations
    n_types = Number of hydrometeor species
    """
    label = ['Drizzle', 'Rain', 'Ice_Crystals', 'Aggregates',
             'Wet_Snow', 'Vertical_Ice', 'Low-Density_Graupel',
             'High-Density_Graupel', 'Hail', 'Big_Drops']
    # 'Large_Hail','SmHail_Rain','LgHail_Rain']
    # DUMMY ARRAYS TO PLOT OVER
    dz = np.arange(701)/10.
    zdr = -2.0 + np.arange(801)/100.
    kdp = -2.0 + np.arange(901)/100.
    ldr = -60.0 + np.arange(5900)/100.
    rho_hv = 0.75 + np.arange(2501)/10000.
    T = -50. + np.arange(1001)/10.
    for i in range(n_types):
        dz_beta = _get_beta_example(dz, mbf_sets, 'Zh_set', i)
        kdp_beta = _get_beta_example(kdp, mbf_sets, 'Kdp_set', i)
        zdr_beta = _get_beta_example(zdr, mbf_sets, 'Zdr_set', i)
        ldr_beta = _get_beta_example(ldr, mbf_sets, 'LDR_set', i)
        rho_beta = _get_beta_example(rho_hv, mbf_sets, 'rho_set', i)
        T_beta = _get_beta_example(T, mbf_sets, 'T_set', i)
        fig, ax = plt.subplots(3, 2, figsize=(8, 5))
        fig.subplots_adjust(top=0.90, bottom=0.1, left=0.05, right=0.95,
                            hspace=0.45)
        ax = ax.flatten()
        fig.suptitle('%s band %s, temp_factor: %.1f' % (band,
                     label[i], temp_factor))
        ax[0].plot(dz, dz_beta, linewidth=2)
        ax[0].set_title('dBZ')
        ax[1].plot(kdp, kdp_beta, linewidth=2)
        ax[1].set_title('Kdp')
        ax[2].plot(zdr, zdr_beta, linewidth=2)
        ax[2].set_title('Zdr')
        ax[3].plot(ldr, ldr_beta, linewidth=2)
        ax[3].set_title('LDR')
        ax[3].set_ylim(0, 1)
        ax[4].plot(rho_hv, rho_beta, linewidth=2)
        ax[4].set_title(r'$\rho$$_{hv}$')
        ax[5].plot(T, T_beta, linewidth=2)
        ax[5].set_title('T')
        plt.savefig(pdir+'%s_beta_%s.png' % (label[i], band))
        plt.close()


def get_beta_set(filename, factor=1.0):
    """
    Function to read the CSVs that contain the beta function parameters for
    each variable and frequency band
    filename = Name of file
    factor = Modifier to effect of variable (normally only used for T)
    """
    beta_parameters = pd.read_csv(filename, index_col=0)
    set_values = {}
    for label in ['m', 'a', 'b']:
        set_values[label] = np.array(beta_parameters.loc[:, label]) / factor
    return set_values

################################


def get_mbf_sets_summer(use_temp=True, plot_flag=False, n_types=10,
                        temp_factor=1, band='S', fdir=CSV_DIR, verbose=False,
                        pdir='./'):
    """
    Define structures to hold the beta function parameters for each input
    variable and fuzzy set. Valid for summer precip.

    HID types:           Fuzzy set:
    -------------------------------
    Drizzle                  1
    Rain                     2
    Ice Crystals             3
    Aggregates               4
    Wet Snow                 5
    Vertical Ice             6
    Low-Density Graupel      7
    High-Density Graupel     8
    Hail			         9
    Big Drops		         10

    Arguments:
    use_temp = Flag to create T beta function sets
    plot_flag = Flag to turn on optional beta function plots
    band = 'C' or 'S'
    temp_factor = Factor to modify depth of T effects; > 1 will broaden the
                  slopes of T MBFs
    fdir = Path to CSV locations
    n_types = Number of hydrometeor species
    verbose = Set to True to get text updates
    pdir = Directory to stash beta function plots (if made)

    Returns:
    mbf_sets = Membership beta function sets
    """

    me = '-->get_mbf_sets_summer  '
    if verbose:
        print(me + ' 10 Category Summer %s-band HID' % band)

    # Horizontal Reflectivity (Zh)
    fname = band + '-band_Reflectivity.csv'
    Zh_set = get_beta_set(fdir+fname)

    # Differential Reflectivity (Zdr)
    fname = band + '-band_Differential_Reflectivity.csv'
    Zdr_set = get_beta_set(fdir+fname)

    # Specific Differential Phase (Kdp)
    fname = band + '-band_Specific_Differential_Phase.csv'
    Kdp_set = get_beta_set(fdir+fname)

    # Linear depolarization ratio (LDR)
    fname = band + '-band_Linear_Depolarization_Ratio.csv'
    LDR_set = get_beta_set(fdir+fname)

    # Correlation coefficient at zero lag (rho)
    fname = band + '-band_Correlation_Coefficient.csv'
    rho_set = get_beta_set(fdir+fname)

    # Temperature (T)
    if use_temp:
        fname = band + '-band_Temperature.csv'
        T_set = get_beta_set(fdir+fname, factor=temp_factor)
        if verbose:
            print(me+'Using Temperature in HID')
    else:
        T_set = None
        if verbose:
            print(me+'No Temperature used in HID')

    # Now populate the full mbf_sets dictionary
    mbf_sets = {'Zh_set': Zh_set, 'Zdr_set': Zdr_set, 'Kdp_set': Kdp_set,
                'LDR_set': LDR_set, 'rho_set': rho_set, 'T_set': T_set}

    if plot_flag:
        plot_beta_functions(mbf_sets, band, temp_factor, pdir=pdir,
                            n_types=n_types)

    return mbf_sets

################################

if __name__ == '__main__':
    a = cdf_betas_sum(plot_flag=True, use_temp=True, temp_factor=2, band='C')
