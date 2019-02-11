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

If using Anaconda/Miniconda, CSU_RadarTools can be installed from `conda-forge` (in progress) with:

```
conda install -c defaults -c conda-forge csu_radartools
```

We suggest you do so in a dedicated conda environment (e.g. shared with PyART and/or your other radar data manipulation tools) to avoid any chance of contaiminating your `base` environment.

Alternatively, if not using Anaconda/Miniconda, you can install CSU_Radartools using `pip` (or `pip3`, if `pip` does not point to your Python 3 install on your system):

```
pip install csu_radartools
```

Again, we suggest using a proper `virtualenv`/`venv` for your radar analysis stack, *particularly* if using your system Python install.
If are using your system Python or one installed for all users, you may need to have sudo privileges, depending on your setup.

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

CSU_RadarTools is known to work under Python 2.7 and 3.4-3.7.
Other Python versions are untested.
Although we still fully support Python 2.7 at this time, we strongly recommend you move to Python >=3.6 as soon as practicable, given the Python 2 end of life date in less than a year, and the fact that many other scientific packages have dropped it already.
For more information, please see e.g. [the Python 3 Statement](https://python3statement.org/).

Latest release of CSU_RadarTools (v1.3):
[![DOI](https://zenodo.org/badge/31606116.svg)](https://zenodo.org/badge/latestdoi/31606116)
