#!/usr/bin/python

from distutils.core import setup
from Cython.Build import cythonize
import numpy as np
setup(
	name='co_var',
	ext_modules=cythonize('./src/Cython/co_var.pyx'),
	include_dirs = [np.get_include()], 
)

setup(
	name='lpa_init',
	ext_modules=cythonize('./src/Cython/lpa_init.pyx'),
	include_dirs = [np.get_include()], 
)
