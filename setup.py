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
import ast
import setuptools

try:
    import numpy
except ImportError:
    raise RuntimeError("Cannot find NumPy. Please install it first "
                       "or use pip >= 10, which will do so automatically.")


# Set CSU_F2PY to True to use f2py instead of Cython for csu_kdp, etc
if ('CSURT_F2PY' in os.environ
        and os.environ.get('CSURT_F2PY').lower() not in {'0', 'false', 'no'}):
    USE_CYTHON = False
else:
    USE_CYTHON = True

if USE_CYTHON:
    from setuptools import setup, Extension  # analysis:ignore
    try:
        import Cython  # analysis:ignore
    except ImportError:
        raise RuntimeError("Cannot find Cython. Please install it first "
                           "or use pip >= 10, which will do so automatically.")
    from Cython.Build import cythonize  # analysis:ignore
    from Cython.Compiler import Options  # analysis:ignore
    Options.language_level = '2'
else:
    from numpy.distutils.core import setup, Extension  # analysis:ignore


# Get current location
HERE = os.path.abspath(os.path.dirname(__file__))
# Pull the header into a variable
DOCLINES = __doc__.split('\n')
# Get packages
PACKAGES = setuptools.find_packages()


# Get package version
def get_version(module=None):
    """Get version string for package."""
    if module is None:
        module = setuptools.find_packages()[0]
    with open(os.path.join(HERE, module, '_version.py'), 'r') as version_file:
        data = version_file.read()
    lines = data.split('\n')
    for line in lines:
        if line.startswith('VERSION_INFO'):
            version_tuple = ast.literal_eval(line.split('=')[-1].strip())
            version = '.'.join(map(str, version_tuple))
            break
    return version


if USE_CYTHON:
    EXT = '.pyx'
else:
    EXT = '.f'

EXTENSIONS = [Extension(PACKAGES[0] + '.calc_kdp_ray_fir',
                        [PACKAGES[0] + '/calc_kdp_ray_fir' + EXT])]
INCLUDE_DIRS = [numpy.get_include(), '.']

if USE_CYTHON:
    EXTENSIONS = cythonize(EXTENSIONS, compiler_directives={'cpow': True})


# Run setup
setup(name='csu_radartools',
      version=get_version(),
      url='https://radarmet.atmos.colostate.edu',
      download_url='https://github.com/CSU-Radarmet/CSU_RadarTools/releases',
      author='Brenda Dolan, Brody Fuchs, Timothy Lang',
      author_email='bdolan@colostate.edu',
      description=DOCLINES[1],
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
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering :: Atmospheric Science',
          ],
      packages=PACKAGES,
      package_data={'csu_radartools': ['beta_function_parameters/*.csv']},
      ext_modules=EXTENSIONS,
      include_dirs=INCLUDE_DIRS,
      install_requires=['numpy', 'pandas', 'matplotlib', 'scipy', 'cython', 'netCDF4'],
      python_requires='<4',
      )
