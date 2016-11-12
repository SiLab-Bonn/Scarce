# Scarce
[![Build Status](https://travis-ci.org/SiLab-Bonn/Scarce.svg?branch=master)](https://travis-ci.org/SiLab-Bonn/Scarce)
[![Coverage Status](https://coveralls.io/repos/github/SiLab-Bonn/Scarce/badge.svg?branch=master)](https://coveralls.io/github/SiLab-Bonn/Scarce?branch=master)

_Scarce_ stands for Silicon ChARge Collection Efficiency and is a software
to calculate the charge collection efficiency of irradiated silicon sensors.
Planar and 3D electrode configuration are supported.
Additionally a large collection of (semi-) empiric formulars is provided to
calculate silicon properties.

# Installation
So far the installation has been tested with the following environment only:
- Ubuntu 14.04 LTS 64-bit
- Anaconda Python 2.7

1. Install the mesh creator gmsh:
`sudo apt-get install gmsh`

2. Install Anaconda Python distribution: 
`wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh` 
For more information visit: https://www.continuum.io/downloads

3. Install precompiled dependencies: 
`conda install numpy pytables scipy matplotlib`

4. Install dependencies:
`pip install -e git://pysparse.git.sourceforge.net/gitroot/pysparse/pysparse#egg=PySparse ez_setup`

5. Download Scarce:
`git checkout https://github.com/SiLab-Bonn/Scarce`

6. Install _Scarce_ in Python develop mode:
`cd Scarce && python setup.py develop`




