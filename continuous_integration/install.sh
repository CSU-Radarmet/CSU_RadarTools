#!/bin/bash
# This script was adpated from the pyart install.sh script.
# This script is adapted from the install.sh script from the scikit-learn
# project: https://github.com/scikit-learn/scikit-learn

# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

set -e
# use next line to debug this script
#set -x

# Use Miniconda to provide a Python environment.  This allows us to perform
# a conda based install of the SciPy stack on multiple versions of Python
# as well as use conda and binstar to install additional modules which are not
# in the default repository.
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh \
    -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b
export PATH=/home/travis/miniconda2/bin:/home/travis/miniconda/bin:$PATH
conda update --yes conda

# Create a testenv with the correct Python version
conda create -n testenv --yes pip python=$TRAVIS_PYTHON_VERSION
source activate testenv

# Install dependencies
conda install --yes pandas cython basemap numpy scipy matplotlib netcdf4 nose hdf4=4.2.12 
pip install skewt pytest
conda install --yes -c conda-forge arm_pyart



#if [[ $PYTHON_VERSION == '2.7' ]]; then
#    pip install xmltodict
#fi

# install coverage modules

pip install -v -e .
