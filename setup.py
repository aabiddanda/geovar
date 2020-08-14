#!/usr/bin/env python

from setuptools import setup

version = '0.1.0'

required = open('requirements.txt').read().split('\n')

setup(
    name='geovar',
    version=version,
    description='A package to provide functions to visualize joint allele frequency spectra[D[D[D[D[D[D[D[D[D[D[D[D[D',
    author='Arjun Biddanda',
    author_email='aabiddanda@gmail.com',
    url='https://github.com/Arjun Biddanda/geovar',
    packages=['geovar'],
    install_requires=required,
    long_description='See ' + 'https://github.com/Arjun Biddanda/geovar',
    license='MIT'
)
