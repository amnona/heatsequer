#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2016--, heatsequer development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re
import ast
import os

from setuptools import find_packages, setup


classes = """
    Development Status :: 3 - Alpha
    License :: OSI Approved :: BSD license
    Topic :: Software Development :: Libraries
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3 :: Only
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

description = ('Easy heatmap visualization of microbes')

with open('README.md') as f:
    long_description = f.read()


version = '0.0.1'

setup(name='heatsequer',
      version=version,
      license='Modified BSD',
      description=description,
      long_description=long_description,
      author="heatsequer development team",
      author_email="amnonim@gmail.com",
      maintainer="heatsequer development team",
      maintainer_email="amnonim@gmail.com",
      packages=find_packages(),
      setup_requires=['numpy >= 1.9.2'],
      # note, hdf5 is required to be installed beforehand
      # also pyqt==4.11.4 needs to be installed beforehand
      install_requires=[
          'biom-format',
          'easygui',
          'scipy',
          'numpy',
          'networkx',
          'scikit-learn',
          'matplotlib',
          'h5py',
          'requests',
      ],
      classifiers=classifiers,
      package_data={
          }
      )
