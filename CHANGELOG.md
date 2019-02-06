# Change Log


## 2016-11

### CSU_Blended_Rain_Tropical

* Added the dBZ + Kdp threshold back in to eliminate small reflectivities
  with non-zero Kdp blowing up the rain rates.



## 2016-07-28

### CSU_DSD

* Changed import location for _check_for_array


### CSU_Blended_Rain

* Moved base rainfall functions to common, and also now import
  _check_for_array from common.



## 2016-05-10

### Common

* Moved functions from csu_blended_rain here.


### CSU_Blended_Rain

* Added ability for user to provide custom parameters to
  polarimetric rainfall equations via the blended rainfall
  routines. Also customized blended routines to handle
  non-S bands. In this case, only R-Z and R-Kdp are used.



## 2016-05-03

### Beta_Functions

* Python 3 compatible.


### CSU_KDP

* Cython now an option for speeding up KDP routines.


### Calc_KDP_Ray_Fir

* Updated for Cython. Confirmed ~50 times faster than 1-deg np.polyfit().


### CSU_FHC

* Cython now an option for speeding up the hid_beta routines.



## 2015-12-08

### CSU_KDP

* Made number of gates used for standard deviation calculation adjustable.



## 2015-11-20


### CSU_KDP

* Now using a Fortran shared object (calc_kdp_ray_fir) to
  do the ray-based KDP calculations. This has vastly sped up
  the overall KDP processing (> 100x). f2py FTW!
* Updated calc_kdp_bringi to fail softly when (window/gs) is not even.


### CSU_FHC

* Sped up hid_beta by using f2py + working w/ 1-D flattened arrays
  that are later reshaped to the necessary shape.


### CSU_Liquid_Ice_Mass

* Performance improvements.



## 2015-09-16

### CSU_Blended_Rain

* Fixed logical inconsistencies leading to lack of
  rainfall calculation in HID = rain + low Z + high Kdp/Zdr.



## 2015-09-04

### CSU_KDP

* Added window keyword to enable stretching the FIR window (e.g.,
  use a 21-pt filter over 5 km with 250-m gate spacing).
* Forcing FIR order to be even, _calc_kdp_ray will crash otherwise.


### CSU_Misc

* Vastly sped up despeckle routine using scipy.



## 2015-08-05

### CSU_KDP

* Made Python 3 compatible.
* Fixed issue with non-integer array indices.


## CSU_FHC

* Python 3.


### CSU_Blended_Rain

* Python 3.


### CSU_Liquid_Ice_Mass

* Python 3 compliant.


### CSU_Misc

* Made Python 3 compatible.
* Made pep8 compatible.



## 2015-07-10

### CSU_KDP

* Made sub-module pep8 compliant.
* Added despeckle() along with a private helper function.
* Added warnings.warn import.



## 2015-05-08

### CSU_Misc

* Added despeckle() along with a private helper function.
* Added warnings.warn import.



## 2015-04-27

### CSU_KDP

* Made algorithm work with a user-defined gate spacing (via gs keyword).
  Untested on gate spacings that do not divide evenly into the 3-km window
  used for filtering the PHIDP data, however. But common gate spacings
  like 50, 100, 150, 200, 250, and 300 meters should all work fine.
* Made the algorithm capable of receiving 2D array inputs (i.e., azimuth &
  range) as well as 1D inputs (range only). If 2D, rng needs to be 2D as
  well. However, thsd should remain a scalar, or 1D and only vary by range.


## Legacy changelog for Beta_Functions:

### 25 January 2015

* Pythonized


### 10 June 2012

* Adjusted MBFs for performance.


### 08 February 2012

* Changed the categories and the MBF values based on scattering simulations.
  NOTE: LDR values from simulations were added.


### 29 January 2009
* Changed the MBFS to the theory-based S-band values.
  NOTE: LDR VALUES WERE NOT MODIFIED.


### 19 February 2003

* Changed MBFs for wet/dry graupel and Wet Snow to
  conform to Liu and Chandrasekar (2000).
* Changed MBF for vertical ice to
  conform to Carey and Rutledge (1998).

### 15 May 2002

* Changed KDP MBF for vertical ice to be unity from -0.2 to -0.6.
  Used to be from -0.6 to 0.0. This didn't make much difference.


### 28 April 2002

* Changed the fuzzy sets to include vertical ice
  and get rid of low/high density snow.







