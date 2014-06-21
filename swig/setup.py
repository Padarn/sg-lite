#!/usr/bin/env python

from distutils.core import setup, Extension

import numpy
import os
import sys

numpyinclude = numpy.__file__[:-12]+'core/include/'
libpath = os.path.realpath('../sparsegrid/')
platform = sys.platform
if platform != 'darwin':
    runtime = '-rpath ' + libpath
else:
    runtime = ''

sparsepad = Extension('_sparsepad',
                           sources=['sparsepad_wrap.cxx'],
                           libraries=['sg-lite'],
                           include_dirs=['../',numpyinclude],
                           library_dirs=[libpath],
                           extra_objects=[runtime]
                           )

setup (name = 'sparsepad',
       version = '0.1',
       author      = "Padarn Wilson",
       description = """sg-lite swig interface""",
       ext_modules = [sparsepad],
       py_modules = ["sparsepad"],
       )