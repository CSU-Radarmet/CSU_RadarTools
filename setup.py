# -*- coding: utf-8 -*-
"""
Python tools for polarimetric radar retrievals.

To access, use the following in your analysis code:
from csu_radartools import (
    csu_fhc, csu_kdp, csu_dsd, csu_liquid_ice_mass, csu_misc,
    csu_blended_rain, fundamentals)

Works on Windows, but you'll need to install a C++ compiler (e.g. MSVC >=2015).
"""

import os
import numpy
import setuptools


# Set to False to use f2py instead of Cython for csu_kdp, etc
if os.environ.get('USE_CYTHON', False):
    USE_CYTHON = False
else:
    USE_CYTHON = True

if USE_CYTHON:
    from setuptools import setup, Extension
    from Cython.Build import cythonize
    from Cython.Compiler import Options
    Options.language_level = '2'
else:
    from numpy.distutils.core import setup, Extension


# Pull the header into a variable
doclines = __doc__.split('\n')

VERSION = '1.2'

# Set variables for setup
PACKAGES = setuptools.find_packages()

if USE_CYTHON:
    EXT = '.pyx'
else:
    EXT = '.f'

extensions = [Extension(PACKAGES[0] + '.calc_kdp_ray_fir',
                        [PACKAGES[0] + '/calc_kdp_ray_fir' + EXT])]

if USE_CYTHON:
    extensions = cythonize(extensions)


# Run setup
setup(name='csu_radartools',
      version=VERSION,
      url='https://radarmet.atmos.colostate.edu',
      download_url='https://github.com/CSU-Radarmet/CSU_RadarTools/releases',
      author='Brenda Dolan, Brody Fuchs, Timothy Lang',
      author_email='bdolan@atmos.colostate.edu',
      description=doclines[1],
      long_description=__doc__,
      keywords='radar precipitation meteorology weather',
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Education',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: POSIX :: Linux',
          'Operating System :: Unix',
          'Programming Language :: Cython',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering :: Atmospheric Science',
          ],
      packages=PACKAGES,
      package_data={'csu_radartools': ['beta_function_parameters/*.csv']},
      ext_modules=extensions,
      include_dirs=[numpy.get_include(), '.'],
      install_requires=['numpy', 'pandas', 'matplotlib', 'scipy', 'cython'],
      python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4',
      )
