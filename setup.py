"""
CSU_RadarTools
"""

import os
import sys

from numpy.distutils.core import setup, Extension

# Pull the header into a variable
doclines = __doc__.split("\n")

VERSION = '1.1'

# Set variables for setup
PACKAGES = ['csu_radartools']

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
          Programming Language :: Python",
          Topic :: Scientific/Engineering
          Topic :: Scientific/Engineering :: Atmospheric Science
          Operating System :: Unix
          Operating System :: POSIX :: Linux
          Operating System :: MacOS
          """],
      long_description="""
          Python tools for polarimetric radar retrievals.
          To access, use the following in your analysis code:
          from csu_radartools import csu_fhc
          """,)
#       ext_modules = [Extension('calc_kdp_ray_fir',
#                      [PACKAGES[0]+'/calc_kdp_ray_fir.f'])],
