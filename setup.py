# -*- coding: utf-8 -*-
"""
 *  File:        setup.py
 *  Author:      Modestino Carbone
 *  Created on:  10/09/2025
 *  Last update: 03/10/2025
 *
 *  Description:
 *      - 'csd_basic'and 'csd_welch' : The main basics of cross-spectral density (CSD).
 *      - 'csd_Min_Max': Cross spectral density (CSD) with
 *      advanced averaging options.
 *
"""
from setuptools import setup, find_packages

setup(
    name="Spectra_lib",
    version="2.1",
    author="Modestino Carbone",
    description="A collection of Python spectrum function for Harmonic analysis",
    packages=find_packages(),
)
