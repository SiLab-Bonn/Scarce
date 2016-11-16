# Scarce
[![Build Status](https://travis-ci.org/SiLab-Bonn/Scarce.svg?branch=master)](https://travis-ci.org/SiLab-Bonn/Scarce)
[![Coverage Status](https://coveralls.io/repos/github/SiLab-Bonn/Scarce/badge.svg?branch=master)](https://coveralls.io/github/SiLab-Bonn/Scarce?branch=master)

_Scarce_ stands for Silicon ChARge Collection Efficiency and is a software
to calculate the charge collection efficiency of irradiated silicon sensors.
Planar and 3D electrode configuration are supported.
Additionally a large collection of (semi-) empiric formulars is provided to
calculate silicon properties.

# Installation

Linux is the preferred operating system, since it supports the fast sparse matrix library pysparse. Nevertheless
one can also get good results in reasonble time under Windows.

## Linux
So far the installation has been tested with the following environment only:
- Ubuntu 14.04 LTS 64-bit with Anaconda Python 2.7

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

6. Install Scarce in development mode by typing: 
`cd Scarce && pip setup.py develop`

## Windows
So far the installation has been tested with the following environment only:
- Windows 7 with Anaconda Python 2.7

1. Install the mesh creator gmsh that can be donwloaded here:
http://gmsh.info/bin/Windows/gmsh-2.14.1-Windows64.zip

2. Install 64-bit Anaconda Python 2.7 distribution that can be donwloaded here:
https://www.continuum.io/downloads#windows

3. Install precompiled dependencies by typing into the command prompt: 
`conda install numpy pytables scipy matplotlib`

4. Download Scarce here and unpack to a folder of your choise:
https://github.com/SiLab-Bonn/Scarce/archive/master.zip

5. Install Scarce in development mode by typing: 
`pip setup.py develop`




