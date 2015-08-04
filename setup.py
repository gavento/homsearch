#!/usr/bin/env python
from distutils.core import setup, Extension
from Cython.Build import cythonize

extension = Extension(
           "homsearch_interface",                 # module name
           sources=["homsearch_interface.pyx", "homsearch_lib.cpp"],  # source files
           language="c++",             # generate C++ code
           extra_compile_args=['-std=c++11'],
           )
    

setup(ext_modules = cythonize(extension))
