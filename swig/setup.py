#!/usr/bin/env python

from distutils.core import setup, Extension

import numpy
import os

numpyinclude = numpy.__file__[:-12]+'core/include/'
#libfullpath = os.path.realpath('../sparsegrid/libsg-lite.so')
libpath = os.path.realpath('../sparsegrid/')


sparsepad = Extension('_sparsepad',
                           sources=['sparsepad_wrap.cxx'],
                           #extra_objects=[libfullpath],
                           libraries=['sg-lite'],
                           include_dirs=['../',numpyinclude],
                           library_dirs=[libpath],
                           )

setup (name = 'sparsepad',
       version = '0.1',
       author      = "Padarn Wilson",
       description = """sg-lite swig interface""",
       ext_modules = [sparsepad],
       py_modules = ["sparsepad"],
       )