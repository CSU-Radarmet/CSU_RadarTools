# CSU_RadarTools

Python tools for polarimetric radar retrievals.

This codebase was developed at Colorado State University by numerous people, including Brenda Dolan, Brody Fuchs, Kyle Wiens, Rob Cifelli, Larry Carey, Timothy Lang, and others.

Currently, fuzzy-logic-based hydrometeor identification, blended rainfall,
DSD retrievals, and liquid/ice mass calculations are supported.
There is also an algorithm that uses a finite impulse response (FIR) filter to process differential phase and calculate specific differential phase.
Finally, there are some tools to do rudimentary QC on the data.

These are supplied as standalone functions that take polarimetric radar data as arguments.
Scalars and arrays are supported as function inputs. The main exception is `csu_kdp.calc_kdp_bringi()` which requires individual rays, sweeps, or volumes of radar data.


## CSU_RadarTools Installation

To install CSU_RadarTools, download the source code and then run the following in the code's home directory:

```pip install .```

By default, if the `CSURT_F2PY` environment variable is not set (or set to a case-insensitive match for {`0`, `false`, `no`}, the package will be compiled using Cython to speed up KDP and other routines.
This enables the widest cross-platform support.
However, if you are able to compile programs using `f2py`, you may want to try setting `CSU_F2PY` to any non-falsy value (e.g. `1`).
This should provide a modest performance improvement.
Under Windows, `f2py` will be difficult to get working, but we've been successful using it under Mac and Linux.

*Note:* If you have previously installed `CSU_Radartools` you may have to completely remove it
from its installation location to get the latest version to work right, since the KDP
routines have been substantially altered.

To easily access all of the tools, you can use the following in your analysis code:

```python
from csu_radartools import (csu_fhc, csu_liquid_ice_mass, csu_blended_rain,
                            csu_dsd, csu_kdp, csu_misc, fundamentals)
```

For help information do help on individual modules.
There is also a demonstration IPython notebook in the notebooks directory that covers all the modules.

CSU_RadarTools is known to work under Python 3.X

Latest release of CSU_RadarTools (v1.4):
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13975472.svg)](https://doi.org/10.5281/zenodo.13975472)

## References for specific modules

### Summer HID:

* Tessendorf, S. A., Miller, L. J., Wiens, K. C., & Rutledge, S. A. (2005). The 29 June 2000 supercell observed during STEPS. Part I: Kinematics and microphysics. Journal of the Atmospheric Sciences, 62(12), 4127-4150.

   (original S-band algorithm)
      
* Dolan, B., & Rutledge, S. A. (2009). A theory-based hydrometeor identification algorithm for X-band polarimetric radars. Journal of Atmospheric and Oceanic Technology, 26(10), 2071-2088.

   (original X-band algorithm)
      
* Dolan, B., Rutledge, S. A., Lim, S., Chandrasekar, V., & Thurai, M. (2013). A robust C-band hydrometeor identification algorithm and application to a long-term polarimetric radar dataset. Journal of Applied Meteorology and Climatology, 52(9), 2162-2186.

   (S-band MBFs and X-band for hail/Ws are derived from the same scattering simulations described within)


### Winter HCA:

* Thompson, E. J., Rutledge, S. A., Dolan, B., Chandrasekar, V., & Cheong, B. L. (2014). A dual-polarization radar hydrometeor classification algorithm for winter precipitation. Journal of Atmospheric and Oceanic Technology, 31(7), 1457-1481.

### Trapezoidal Temperature functions:

* Rutledge, Steven A., V. Chandrasekar, Brody Fuchs, Jim George, Francesc Junyent, Brenda Dolan, Patrick C. Kennedy, and Kyla Drushka. "SEA-POL Goes to Sea." Bulletin of the American Meteorological Society 100, no. 11 (2019): 2285-2301.

### Kdp calculation

* Hubbert, J., and V. N. Bringi. "An iterative filtering technique for the analysis of copolar differential phase and dual-frequency radar measurements." Journal of Atmospheric and Oceanic Technology 12.3 (1995): 643-648. 


### Ice and liquid mass calculation

* Cifelli, R., Petersen, W. A., Carey, L. D., Rutledge, S. A., & da Silva Dias, M. A. (2002). Radar observations of the kinematic, microphysical, and precipitation characteristics of two MCSs in TRMM LBA. Journal of Geophysical Research: Atmospheres, 107(D20), LBA-44.


### CSU Blended Rainfall Algorithm

* Cifelli, R., Chandrasekar, V., Lim, S., Kennedy, P. C., Wang, Y., & Rutledge, S. A. (2011). A new dual-polarization radar rainfall algorithm: Application in Colorado precipitation events. Journal of Atmospheric and Oceanic Technology, 28(3), 352-364.


### CSU Tropical Rainfall Algorithm

* Thompson, E. J., Rutledge, S. A., Dolan, B., Thurai, M., & Chandrasekar, V. (2018). Dual-polarization radar rainfall estimation over tropical oceans. Journal of Applied Meteorology and Climatology, 57(3), 755-775.



### DSD calculations

* Bringi, V. N., Huang, G. J., Chandrasekar, V., & Gorgucci, E. (2002). A methodology for estimating the parameters of a gamma raindrop size distribution model from polarimetric radar data: Application to a squall-line event from the TRMM/Brazil campaign. Journal of Atmospheric and Oceanic Technology, 19(5), 633-645.

   (Alternate S-band retrieval)

* Bringi, V. N., Williams, C. R., Thurai, M., & May, P. T. (2009). Using dual-polarized radar and dual-frequency profiler for DSD characterization: A case study from Darwin, Australia. Journal of Atmospheric and Oceanic Technology, 26(10), 2107-2122.

  (C-band)

* Bringi, V. N., L. Tolstoy, M. Thurai, and W. A. Petersen. "Estimation of spatial correlation of rain drop size distribution parameters and rain rates using NASA’s S-band polarimetric radar and 2D video disdrometer network: Two case studies from MC3E." In 36th Conf. on Radar Meteorology. 2013.

   (S-band retrieval)

