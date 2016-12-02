===============================================
Introduction
===============================================
|travis-status|  |appveyor-status|  |rtd-status|  |coverage|

Scarce stands for **s**\ ilicon **c**\ h\ **ar**\ ge **c**\ ollection **e**\ fficiency and is a software
to calculate the charge collection-efficiency of irradiated and segmented silicon 
sensors. Planar and 3D electrode configurations are supported.
Additionally a collection of formulars is provided to
calculate silicon properties.

Installation
============
The installation works with Linux and Windows. Mac OS might also work.

Linux
-----
This installation has been tested with Ubuntu 14.04 LTS 64-bit and Anaconda Python 2.7 64-bit.

1. Install the mesh creator gmsh:

.. code-block:: bash
   
   sudo apt-get install gmsh

2. Install Anaconda Python distribution:

.. code-block:: bash

   wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh

For more information visit: https://www.continuum.io/downloads

3. Install precompiled dependencies:

.. code-block:: bash

   conda install numpy pytables scipy matplotlib

4. Install sparse matrix solver (optional, increases speed):

.. code-block:: bash

   pip install -e git://pysparse.git.sourceforge.net/gitroot/pysparse/pysparse#egg=PySparse ez_setup

5. Download Scarce:

.. code-block:: bash

   git checkout https://github.com/SiLab-Bonn/Scarce

6. Install Scarce in development mode by typing: 

.. code-block:: bash

   cd Scarce && python setup.py develop

Windows
-------
This installation has been tested with Windows 7 64-bit and Anaconda Python 2.7 64-bit.

1. Install the mesh creator gmsh that can be donwloaded here:

.. code-block:: bash

   http://gmsh.info/bin/Windows/gmsh-2.14.1-Windows64.zip

2. Install 64-bit Anaconda Python 2.7 distribution that can be donwloaded here:
https://www.continuum.io/downloads#windows

3. Install precompiled dependencies by typing into the command prompt:

.. code-block:: bash

   conda install numpy pytables scipy matplotlib

4. Download Scarce here and unpack to a folder of your choise:
https://github.com/SiLab-Bonn/Scarce/archive/master.zip

5. Install Scarce in development mode by typing:

.. code-block:: bash

   python setup.py develop

.. |travis-status| image:: https://travis-ci.org/SiLab-Bonn/Scarce.svg?branch=master
    :target: https://travis-ci.org/SiLab-Bonn/Scarce
    :alt: Build status
    
.. |appveyor-status| image:: https://ci.appveyor.com/api/projects/status/32o1x5kcss45m35d?svg=true
    :target: https://ci.appveyor.com/project/DavidLP/scarce
    :alt: Build status

.. |rtd-status| image:: https://readthedocs.org/projects/scarce/badge/?version=latest
    :target: http://scarce.rtfd.org
    :alt: Documentation
    
.. |coverage| image:: https://coveralls.io/repos/github/SiLab-Bonn/Scarce/badge.svg?branch=master
    :target: https://coveralls.io/github/SiLab-Bonn/Scarce?branch=master
    :alt: Coverage

