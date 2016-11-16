#!/usr/bin/env python
from setuptools import setup, find_packages  # This setup relies on setuptools since distutils is insufficient and badly hacked code

version = '0.0.1'
author = 'David-Leon Pohl'
author_email = 'pohl@physik.uni-bonn.de'

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='Scarce',
    version=version,
    description='Charge collection efficiency simulation for irradiated silicon detectors.',
    url='https://github.com/SiLab-Bonn/Scarce',
    license='MIT License',
    long_description='A FVM/FDM based software to calculate the charge collection efficiency of irradiated silicon sensors. Planar and 3D electrode configuration is supported',
    author=author,
    maintainer=author,
    author_email=author_email,
    maintainer_email=author_email,
    packages=find_packages(),
    setup_requires=['setuptools', 'ez_setup'],
    install_requires=required,
    include_package_data=True,  # accept all data files and directories matched by MANIFEST.in or found in source control
    package_data={'': ['README.*', 'VERSION'], 'docs': ['*'], 'examples': ['*']},
    keywords=['silicon', 'detector', 'CCE', 'charge', 'efficiency'],
    platforms='any'
)
