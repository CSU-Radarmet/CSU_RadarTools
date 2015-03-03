"""
csu_fhc.py

Brody Fuchs, CSU, Sept 2014
brfuchs@atmos.colostate.edu

Porting over Brenda Dolan's HID code from IDL 
(apparently originally from Kyle Wiens)

Modifications by Timothy Lang
tjlangoc@gmail.com
1/21/2015

#### BELOW ARE COPIED COMMENTS AND INFO FROM ORIGINAL IDL/WIENS CODE #####
#### MANY OF THESE ARE NO LONGER APPLICABLE; PROVIDED ONLY FOR POSTERITY #####
#### UPDATED COMMENTS ARE IN DOC STRINGS FOR PYTHON FUNCTIONS #####
#+
# NAME:
# cdf_fhc
#
# PURPOSE:
#  Perform the fuzzy logic hydrometeor classification (FHC)
#  using input polarimetric variables and temperature.  Return the FHC
#  types in an array of the same size as the input polarimetric variables.
#
#  This method is based on the method of Liu and Chandrasekar (2000),
#  J. Atmos. Oceanic Tech., _17_, 140-164.  This assumes summer convective
#  storms.   
#  However, my method uses a weighted sum of each truth value instead
#  of the straight multiplication method of L&C (2000).
#
# Truth value of each hydrotype for each input variable is a beta
# function.  The total truth value for each hydro type is thus a
# weighted sum of the truth value of each input variable, with two
# exceptions:  The horizontal reflectivity and temperature (if used) truth
# values are multiplied by the weighted sum of the others, thus effectively
# weighting the final output more strongly on the horizontal reflectivity
# and temperature (if present).  
#
# Uses parameters defined in 'cdf_beta_functions.pro' to define the
# m,a,b parameters for each hydro type's beta function.  The independent
# variable 'x' is one of the radar measurements.
#
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
# 
# Anticipate having up to 6 input data fields and 11 fuzzy sets (HID types) for
# summer storms.
#
# Input data fields:           HID types:           Fuzzy set:
#-------------------           ----------           ----------
#     Z_h                      Drizzle                  1
#     Z_dr                     Rain                     2
#     K_dp                     Dry snow                 3
#     LDR                      Wet snow                 4
#     rho_hv[0]                Vertical Ice             5
#     Temperature              Dry graupel              6
#                              Wet graupel              7
#                              Small hail               8
#                              Large hail               9
#                              Sm. Hail & Rain         10
#                              Lg. Hail & Rain         11
#
#  So there should be 6 X 11 = 66 beta functions...one for each pair
#  of inputs and fuzzy sets
#
# The final output for each grid point is the fuzzy set which has the
# highest cumulative truth value (i.e., the most likely hydrometeor type).
#
# If there is an error, the function will return the integer zero instead
# of a gridded array.
#
# USAGE:
# output = cdf_fhc(radar_data,mbf_sets,fhc_vars[,T,weights=weights,
#                  use_temp=use_temp,compare=compare])
#
# Terms in brackets [] are optional.
#
# INPUT PARAMETERS:
#        REQUIRED INPUTS
# radar_data -- A structure containing the radar data 
#               (Should not matter if in radial
#                or gridded format).  The structure must
#               have one of the following tag names: DZ,DR,KD,LD,RH. 
#               These correspond to horiz. reflect, diff. reflect.,
#               spec. diff. phase, linear depol. ratio, and corr. coeff.,
#               respectively.  For example, to use all these variables in
#               the classification, define this structure as:
#
#               radar_data = {DZ:dz_data,$
#                             DR:DR_data,$
#                             KD:KD_data,$
#                             LD:LD_data,$
#                             RH:RH_data}
#
#               Then call the program:
#               result = cdf_fhc(radar_data,...)
#
#               Or if you have only horiz. reflect. and
#               diff. reflect. define the structure like this:
#              
#               radar_data = {DZ:dz_data,$
#                             DR:DR_data}
#
#               Only the tag names are important, e.g., the two-character
#               string before each colon must be DZ, DR, etc.  The
#               actual data (e.g., dz_data) should be 1,2,or3-d floating
#               point array of the appropriate data.  Each data field is
#               expected to have the same dimensions.
# mbf_sets -- Beta function parameters for each fuzzy set and input
#             variable (defined in 'cdf_beta_functions.pro')
# fhc_vars -- A Boolean structure of which input variables to use in the
#             FHC.  Should be a 6-element structure of either 0 (Don't
#             Include) or 1 (Do Include)
#             fhc_vars = {Zh: 0 or 1,
#                         Zdr:0 or 1,
#                         Kdp:0 or 1,
#                         LDR:0 or 1,
#                         rho_hv: 0 or 1,
#                         T: 0 or 1}
#             The use for this array is to do sensitivity tests where you'd
#             explicity state which radar fields to use in the
#             classification.  If the corresponding data is not present
#             in the radar_data input structure, then this program will
#             override what you passed it as 'fhc_vars' and set
#             its boolean value to 0.  For example, if there is no LD field
#             in the radar_data structure, then this program will do the
#             following:
#
#             fhc_vars.LDR = 0
#
#             And hence will not try to use LDR in the classification.
#
#       OPTIONAL INPUTS
# T -- Temperature, either single value or a 1-D vertical profile with the
#      temperature defined at each vertical grid point of the radar data.
#        --BD 2-2012: Changed to require this be the same size as the input data.
# weights -- Weights for the different variables. Added this as a keyword such 
#           that you don't need a speciifc cdf_fhc for different radars.  
#           If weights is not passed, then default values are used. 
#           This is a structure like fhc_vars:
#        Weights= { W_Zh: 1.5 ,
#                   W_Zdr: 0.8,
#                   W_Kdp: 1.0,
#                   W_rho: 0.8,
#                   W_ldr: 0.5,
#                   W_T: 0.4}
#
#
# KEYWORDS:
# use_temp -- If set, use temperature in the FHC.  If you set this keyword,
#             make sure you also supply a temperature (T) input.
# compare -- If set, a structure will be returned intstead of just the FHT
#            variable.  It will contain the first
#             and second best scoring categories, as well as their scores.
#             fdat = { FHT, FHT2, mu_best, mu_sec}
#
#BD 2-2012: CHANGED SO NONE OF THESE ARE AVAILABLE ANYMORE. By making temperature
#            The same size as the other variables, no longer need /rhi or /ppi.
#            Sum_it never quite worked right for me anyway, so stick with the 
#           Hybrid method for now.
# rhi -- IF set, assume temperature is an array of temperatures at
#        each z level of the radar field arrays.  Program assumes that the
#        second dimension of the input radar data corresponds to height.
# ppi -- If set, assume temperature is a single number for a given
#        height level.
# sum_it -- If set, use a weighted sum method for all variables.  If not
#           set, use a weighted sum for Zdr, Kdp, LDR, and rho_hv; then
#           multiply this sum by the scores for Zh and Temperature.
#
# You need to set either the 'rhi' or 'ppi' keyword only if you are passing
# a 2-D cross-section of data to this program and you are including
# temperature.  Otherwise, if you are passing a 3-D array and including
# temperature, the program will assume the 3rd dimension of the input radar
# data corresponds to height.   
#
#IF radar_data is just single values, then the scores for each category given 
#   those inputs will be printed out.
#
#IF the compare keyword is set, then the output will be a structure that contains
#   the first and second best scoring
#    HID types, as well as the socres associated with those two categories.
#
# AUTHOR:
# Kyle C. Wiens
# kwiens@atmos.colostate.edu
#
# DATE:
# 26 April 2002.
#
# MODIFICATIONS:
#
# 28 April 2002:  Changed the aggregation operations.  Instead of
# taking a sum of all the MBFs, I take a sum of the polarimetric MBFs
# and then multiply this sum by the MBFs for temperature and
# horizontal reflectivity.
#
# 4 May 2002:  Added weighting factors.
#              Also added ability to use any combination of the input
#              variables in the FHC by multiplying each MBF by its
#              associated fhc_vars value.  This is for doing sensitivity
#              tests or if you don't have all 5 of the radar measurements.
#
# 5 May 2002:  Made weighting factors for Zdr, Kdp, LDR, and rho_hv
# dependent on the hydro type.
#
# 15 May 2002: Added division by the sum of the weighting factors so
# that the scores are normalized with respect to these factors.  I
# should have done this a long time ago!!!
#
# 24 May 2002: Changed the weighting factors back to being the same for each
# hydrometeor type.
#
# 23 July 2002: Greatly improved the efficiency of this program by
# performing all the beta function calculations on the entire 2 or 3
# dimensional data sets, all at once.  I used to have 4 embedded FOR
# loops to do this.  It took some brain power to figure out how to
# keep the maximum scores from each iteration, but apparently, my
# brain was up to the task.  This runs MUCH faster now, with identical
# results. 
#
# 24 July 2002: Decomposed rain/hail mix into large and small hail
#               components.
#
# 7 April 2003: Added a few lines to allow the user to do the
#               classification using any combination of input variables.
#
# 23 April 2003:  Changed the calling sequence so that the input radar data
# is a structure instead of each individual field.  This should make the
# program more modular by allowing it to work with any combination of radar
# variables.   I also added more comments.
#
# 21 January 2004:  Added option to use a weighted sum for all variables by
#                   passing the /Sum_It keyword.
# 22 February 2012: Modified to be more modular, and no longer require so many 
#                   keywords.
#                    Also set up to pass in a single set of radar observations to 
#                   calculate the scores.
#                    B.D. 2-2012
#
#--------
# This is a sub-function of the main program.
#
# FUNCTION hid_beta
# Compute the value of the beta function from inputs x,a,b,m
# Where,
#        x = independent variable (DZ,KDP,etc.)
#        a = width.
#        b = slope.
#        m = midpoint
#--------

########## END OF ORIGINAL COMMENTS AND INFORMATION ##########

"""

import numpy as np
from beta_functions import get_mbf_sets_summer

DEFAULT_WEIGHTS = {'DZ': 1.5, 'DR': 0.8, 'KD': 1.0, 'RH': 0.8, 'LD': 0.5, 'T': 0.4}

def hid_beta(x_arr, a, b, m):
    """Beta function calculator"""
    return 1.0/(1.0 + ( ((x_arr - m)/a)**2)**b)

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
    band = 'C' or 'S'
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
    mu = Input array + additional dimension containing weights for each HID species
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
        print 'FHC fail, no reflectivity field'
        return None
    if T is None:
        use_temp = False

    #Populate fhc_vars and radar_data based on what was passed to function
    radar_data, fhc_vars = _populate_vars(dz, zdr, kdp, rho, ldr, T, verbose)

    #Now grab the membership beta function parameters
    mbf_sets = get_mbf_sets_summer(use_temp=use_temp, plot_flag=plot_flag,
               n_types=n_types, temp_factor=temp_factor, band=band,
               verbose=verbose)
    sets = _convert_mbf_sets(mbf_sets)

    #Check for presence of polarimetric variables
    pol_flag = _get_pol_flag(fhc_vars)

    #Check for presence of temperature
    if use_temp:
        if verbose:
            print 'Using T in FHC'
    else:
        fhc_vars['T'] = 0
        if verbose:
            print 'Not using T in FHC'

    #Get weighted sums
    weight_sum, varlist = _get_weight_sum(fhc_vars, weights, method, verbose)
    if weight_sum is None:
        return None

    #Now loop over every hydrometeor class
    test_list = _get_test_list(fhc_vars, weights, radar_data, sets, varlist,
                               weight_sum, pol_flag, use_temp, method)
    if test_list is None:
        return None

    #Finish up
    mu = np.array(test_list)
    if verbose:
        print mu.shape
        print 'mu max: ', mu.max()
    return mu

#########################
#Private Functions Below#
#########################

def _convert_mbf_sets(mbf_sets):
    """Simply gets mbf_sets dict into form that matches labels used in csu_fhc"""
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
    """Check for presence of each var, and update dicts as needed"""
    fhc_vars = {'DZ': 1, 'DR': 1, 'KD': 1, 'LD': 1, 'RH': 1, 'T': 1}
    radar_data = {}
    if dz is not None:
        radar_data['DZ'] = dz
    if zdr is not None:
        radar_data['DR'] = zdr
    if kdp is not None:
        radar_data['KD'] = kdp
    if rho is not None:
        radar_data['RH'] = rho
    if ldr is not None:
        radar_data['LD'] = ldr
    if T is not None:
        radar_data['T'] = T
    if 'DR' not in radar_data.keys():
        fhc_vars['DR'] = 0
    if 'KD' not in radar_data.keys():
        fhc_vars['KD'] = 0
    if 'LD' not in radar_data.keys():
        fhc_vars['LD'] = 0
    if 'RH' not in radar_data.keys():
        fhc_vars['RH'] = 0
    if 'T' not in radar_data.keys():
        fhc_vars['T'] = 0
    if verbose:
        print 'USING VARIABLES: ', fhc_vars
    return radar_data, fhc_vars

def _get_weight_sum(fhc_vars, weights, method, verbose):
    """Gets sum of weights and varlist used, which depend on method"""
    if 'hybrid' in method:
        if verbose:
            print 'Using hybrid HID method. Pol vars weighted,',\
                      'Z and T (if used) are multiplied'
        varlist = ['DR', 'KD', 'RH', 'LD']
    elif 'linear' in method:
        if verbose:
            print 'NOT using hybrid, all variables treated as weighted sum'
        varlist = ['DR', 'KD', 'RH', 'LD', 'T', 'DZ']
    else:
        print 'No weighting method defined, use hybrid or linear'
        return None, None
    weight_sum = np.sum(np.array([fhc_vars[key]*weights[key]
                                 for key in varlist]))
    if verbose:
        print 'weight_sum: ', weight_sum
    return weight_sum, varlist

def _calculate_test(fhc_vars, weights, radar_data, sets, varlist, weight_sum, c):
    """Loop over every var to get initial value for each HID species 'test'"""
    test = (np.sum(np.array([fhc_vars[key]*weights[key]*\
            hid_beta(radar_data[key], sets[key]['a'][c],\
            sets[key]['b'][c], sets[key]['m'][c])\
            for key in varlist if key in radar_data.keys()]),\
            axis = 0))/weight_sum
    return test

def _get_test_list(fhc_vars, weights, radar_data, sets, varlist, weight_sum,
                   pol_flag, use_temp, method):
    """
    Master loop to compute HID values for each species ('test' & 'test_list').
    Depending on method used, approach is modfied.
    Currently disabling testing as it gets spoofed by bad data. Better to let the
    calculations continue then mask out the bad data using other methods.
    TO DO: Change poor naming scheme for variables 'test' and 'test_list']
    """
    #TJL - Check order of if statements
    test_list = []
    for c in range(len(sets['DZ']['m'])):
        if 'hybrid' in method: #Hybrid emphasizes Z and T extra HARD
            if pol_flag:
                test = _calculate_test(fhc_vars, weights, radar_data, sets,
                                       varlist, weight_sum, c)
                #if test.max() > 1: #Maximum of test should never be more than 1
                #    print 'Fail loc 1, test.max() =', test.max()
                #    return None
            if use_temp:
                if pol_flag:
                    #*= multiplies by new value and stores in test
                    test *= hid_beta(radar_data['T'], sets['T']['a'][c],
                                     sets['T']['b'][c], sets['T']['m'][c])
                    #print 'in loc 2'
                    #if test.max() > 1: #Maximum of test should never be > 1
                    #    print 'Fail loc 2, test.max() =', test.max()
                    #    return None
                else:
                    test = hid_beta(radar_data['T'], sets['T']['a'][c],
                                    sets['T']['b'][c], sets['T']['m'][c])
            if fhc_vars['DZ']:
                if pol_flag or use_temp:
                    test *= hid_beta(radar_data['DZ'], sets['DZ']['a'][c],
                                     sets['DZ']['b'][c], sets['DZ']['m'][c])
                    #if test.max() > 1: #Maximum of test should never be > 1
                    #    print 'Fail loc 3, test.max() =', test.max()
                    #    return None
                else:
                    test = hid_beta(radar_data['DZ'], sets['DZ']['a'][c],
                                    sets['DZ']['b'][c], sets['DZ']['m'][c])
        elif 'linear' in method: #Just a giant weighted sum
            if pol_flag:
                test = _calculate_test(fhc_vars, weights, radar_data, sets,
                                       varlist, weight_sum, c)
        test_list.append(test)
    return test_list
