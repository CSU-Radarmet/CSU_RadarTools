"""
csu_fhc.py

Brody Fuchs, CSU, Sept 2014
brfuchs@atmos.colostate.edu

Porting over Brenda Dolan's HID code from IDL
(apparently originally from Kyle Wiens)

Modifications by Timothy Lang
tjlangoc@gmail.com
01/21/2015
08/05/2015 - Python 3
11/20/2015 - Sped up hid_beta by using f2py + working w/ 1-D flattened arrays
             that are later reshaped to the necessary shape.
05/03/2016 - Cython now an option for speeding up the hid_beta routines.

"""

from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from .beta_functions import get_mbf_sets_summer
from calc_kdp_ray_fir import hid_beta_f

DEFAULT_WEIGHTS = {'DZ': 1.5, 'DR': 0.8, 'KD': 1.0, 'RH': 0.8, 'LD': 0.5,
                   'T': 0.4}


def hid_beta(x_arr, a, b, m):
    """Beta function calculator"""
    return 1.0/(1.0 + (((x_arr - m)/a)**2)**b)


def csu_fhc_summer(use_temp=True, weights=DEFAULT_WEIGHTS, method='hybrid',
                   dz=None, zdr=None, ldr=None, kdp=None, rho=None, T=None,
                   verbose=False, plot_flag=False, n_types=10, temp_factor=1,
                   band='S'):
    """
    Does FHC for warm-season precip.

    Arguments:
    use_temp = Set to False to not use T in HID
    weights = Dict that contains relative weights for every variable; see
              DEFAULT_WEIGHTS for expected stucture
    method = Currently support 'hybrid' or 'linear' methods; hybrid preferred
    verbose = Set to True to get text updates
    plot_flag = Flag to turn on optional beta function plots
    band = 'X', 'C', or 'S'
    temp_factor = Factor to modify depth of T effects; > 1 will broaden the
                  slopes of T MBFs
    n_types = Number of hydrometeor species
    verbose = Set to True to get text updates

    Input measurands (if not None, all must match in shape/size):
    dz = Input reflectivity scalar/array
    zdr = Input reflectivity scalar/array
    ldr = Input reflectivity scalar/array
    kdp = Input reflectivity scalar/array
    rho = Input reflectivity scalar/array
    T = Input temperature scalar/array

    Returns:
    mu = Input array + addtl dimension containing weights for each HID species
    To get dominant species number: fh = np.argmax(mu, axis=0) + 1

    HID types:           Species #:
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

    """

    if dz is None:
        print('FHC fail, no reflectivity field')
        return None
    if T is None:
        use_temp = False

    # Populate fhc_vars and radar_data based on what was passed to function
    radar_data, fhc_vars, shp, sz = \
        _populate_vars(dz, zdr, kdp, rho, ldr, T, verbose)

    # Now grab the membership beta function parameters
    mbf_sets = get_mbf_sets_summer(
        use_temp=use_temp, plot_flag=plot_flag, n_types=n_types,
        temp_factor=temp_factor, band=band, verbose=verbose)
    sets = _convert_mbf_sets(mbf_sets)

    # Check for presence of polarimetric variables
    pol_flag = _get_pol_flag(fhc_vars)

    # Check for presence of temperature
    if use_temp:
        if verbose:
            print('Using T in FHC')
    else:
        fhc_vars['T'] = 0
        if verbose:
            print('Not using T in FHC')

    # Get weighted sums
    weight_sum, varlist = _get_weight_sum(fhc_vars, weights, method, verbose)
    if weight_sum is None:
        return None

    # Now loop over every hydrometeor class
    test_list = _get_test_list(fhc_vars, weights, radar_data, sets, varlist,
                               weight_sum, pol_flag, use_temp, method, sz)
    if test_list is None:
        return None

    # Finish up
    mu = np.array(test_list)
    shp = np.concatenate([[n_types], shp])
    if verbose:
        print(mu.shape)
        print('mu max: ', mu.max())
    # return mu but make sure the shape is an int array
    return mu.reshape(shp.astype(np.int32))

##########################
# Private Functions Below#
##########################


def _convert_mbf_sets(mbf_sets):
    """Gets mbf_sets dict into form that matches labels used in csu_fhc"""
    sets = {}
    sets['DZ'] = mbf_sets['Zh_set']
    sets['DR'] = mbf_sets['Zdr_set']
    sets['KD'] = mbf_sets['Kdp_set']
    sets['LD'] = mbf_sets['LDR_set']
    sets['RH'] = mbf_sets['rho_set']
    sets['T'] = mbf_sets['T_set']
    return sets


def _get_pol_flag(fhc_vars):
    """Check for presence of polarimetric variables"""
    if fhc_vars['DR'] or fhc_vars['KD'] or fhc_vars['LD'] or fhc_vars['RH']:
        pol_flag = True
    else:
        pol_flag = False
    return pol_flag


def _populate_vars(dz, zdr, kdp, rho, ldr, T, verbose):

    """
    Check for presence of each var, and update dicts as needed.
    Flattens multi-dimensional arrays to optimize processing.
    The output array from csu_fhc_summer() will be re-dimensionalized later.
    """
    varlist = [dz, zdr, kdp, rho, ldr, T]
    keylist = ['DZ', 'DR', 'KD', 'RH', 'LD', 'T']
    fhc_vars = {}
    radar_data = {}
    for i, key in enumerate(keylist):
        var = varlist[i]
        if var is not None:
            if key == 'DZ':
                shp = np.shape(var)
                sz = np.size(var)
            if np.ndim(var) > 1:
                radar_data[key] = np.array(var).ravel().astype('float32')
            elif np.ndim(var) == 1:
                radar_data[key] = np.array(var).astype('float32')
            else:
                radar_data[key] = np.array([var]).astype('float32')
            fhc_vars[key] = 1
        else:
            fhc_vars[key] = 0
    if verbose:
        print('USING VARIABLES: ', fhc_vars)
    return radar_data, fhc_vars, shp, sz


def _get_weight_sum(fhc_vars, weights, method, verbose):
    """Gets sum of weights and varlist used, which depend on method"""
    if 'hybrid' in method:
        if verbose:
            print('Using hybrid HID method. Pol vars weighted,',
                  'Z and T (if used) are multiplied')
        varlist = ['DR', 'KD', 'RH', 'LD']
    elif 'linear' in method:
        if verbose:
            print('NOT using hybrid, all variables treated as weighted sum')
        varlist = ['DR', 'KD', 'RH', 'LD', 'T', 'DZ']
    else:
        print('No weighting method defined, use hybrid or linear')
        return None, None
    weight_sum = np.sum(np.array([fhc_vars[key]*weights[key]
                                 for key in varlist]))
    if verbose:
        print('weight_sum: ', weight_sum)
    return weight_sum, varlist


def _calculate_test(fhc_vars, weights, radar_data, sets,
                    varlist, weight_sum, c, sz):
    """Loop over every var to get initial value for each HID species 'test'"""
#    test = (np.sum(np.array([fhc_vars[key] * weights[key] *
#                            hid_beta(radar_data[key], sets[key]['a'][c],
#                            sets[key]['b'][c], sets[key]['m'][c])
    test = (np.sum(np.array([fhc_vars[key] * weights[key] *
                            hid_beta_f(sz, radar_data[key], sets[key]['a'][c],
                            sets[key]['b'][c], sets[key]['m'][c])
            for key in varlist if key in radar_data.keys()]),
            axis=0))/weight_sum
    return test


def _get_test_list(fhc_vars, weights, radar_data, sets, varlist, weight_sum,
                   pol_flag, use_temp, method, sz):
    """
    Master loop to compute HID values for each species ('test' & 'test_list').
    Depending on method used, approach is modfied.
    Currently disabling testing as it gets spoofed by bad data. Letting the
    calculations continue then mask out the bad data using other methods.
    TO DO: Change poor naming scheme for variables 'test' and 'test_list']
    """
    # TJL - Check order of if statements
    test_list = []
    for c in range(len(sets['DZ']['m'])):
        if 'hybrid' in method:  # Hybrid emphasizes Z and T extra HARD
            if pol_flag:
                test = _calculate_test(fhc_vars, weights, radar_data, sets,
                                       varlist, weight_sum, c, sz)
                # if test.max() > 1:  # Max of test should never be > 1
                #     print 'Fail loc 1, test.max() =', test.max()
                #     return None
            if use_temp:
                if pol_flag:
                    # *= multiplies by new value and stores in test
                    test *= hid_beta_f(sz, radar_data['T'], sets['T']['a'][c],
                                       sets['T']['b'][c], sets['T']['m'][c])
                    # print 'in loc 2'
                    # if test.max() > 1: #Maximum of test should never be > 1
                    #     print 'Fail loc 2, test.max() =', test.max()
                    #     return None
                else:
                    test = hid_beta_f(sz, radar_data['T'], sets['T']['a'][c],
                                      sets['T']['b'][c], sets['T']['m'][c])
            if fhc_vars['DZ']:
                if pol_flag or use_temp:
                    test *= hid_beta_f(
                        sz, radar_data['DZ'], sets['DZ']['a'][c],
                        sets['DZ']['b'][c], sets['DZ']['m'][c])
                    # if test.max() > 1:  # Max of test should never be > 1
                    #     print 'Fail loc 3, test.max() =', test.max()
                    #     return None
                else:
                    test = hid_beta_f(sz, radar_data['DZ'], sets['DZ']['a'][c],
                                      sets['DZ']['b'][c], sets['DZ']['m'][c])
        elif 'linear' in method:  # Just a giant weighted sum
            if pol_flag:
                test = _calculate_test(fhc_vars, weights, radar_data, sets,
                                       varlist, weight_sum, c, sz)
        test_list.append(test)
    return test_list
