"""
CSU_RadarTools
"""

import os
import sys
import numpy

# Set to False to use f2py instead of Cython for csu_kdp, etc.
USE_CYTHON = True

if USE_CYTHON:
    from distutils.core import setup
    from distutils.extension import Extension
    from Cython.Build import cythonize
else:
    from numpy.distutils.core import setup, Extension

# Pull the header into a variable
doclines = __doc__.split("\n")

VERSION = '1.2'

# Set variables for setup
PACKAGES = ['csu_radartools']

if USE_CYTHON:
    ext = '.pyx'
else:
    ext = '.f'

extensions = [Extension('calc_kdp_ray_fir',
              [PACKAGES[0]+'/calc_kdp_ray_fir'+ext])]

if USE_CYTHON:
    extensions = cythonize(extensions)

# Run setup
setup(
      name='csu_radartools',
      version=VERSION,
      url='http://radarmet.atmos.colostate.edu',
      author='Brenda Dolan, Brody Fuchs, Timothy Lang',
      author_email='bdolan@atmos.colostate.edu',
      description=doclines[0],
      license='LICENSE',
      packages=PACKAGES,
      package_data={'csu_radartools': ['beta_function_parameters/*.csv']},
      classifiers=["""
          Development Status :: Beta,
          Programming Language :: Python,
          Topic :: Scientific/Engineering
          Topic :: Scientific/Engineering :: Atmospheric Science
          Operating System :: Unix
          Operating System :: POSIX :: Linux
          Operating System :: MacOS
          """],
      long_description="""
          Python tools for polarimetric radar retrievals.
          To access, use the following in your analysis code:
          from csu_radartools import (
              csu_fhc, csu_kdp, csu_dsd, csu_liquid_ice_mass, csu_misc,
              csu_blended_rain, fundamentals)
          """,
      ext_modules = extensions,
      include_dirs = [numpy.get_include()]
      )
