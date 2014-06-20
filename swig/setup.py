#!/usr/bin/env python

from distutils.core import setup, Extension

import numpy
numpyinclude = numpy.__file__[:-12]+'core/include/'

sparsepad = Extension('_sparsepad',
                           sources=['sparsepad_wrap.cxx'],
                           extra_objects=['../sparsegrid/libsg-lite.so'],
                           include_dirs=['../',numpyinclude],
                           library_dirs=['../sparsegrid']
                           )

setup (name = 'sparsepad',
       version = '0.1',
       author      = "Padarn Wilson",
       description = """sg-lite swig interface""",
       ext_modules = [sparsepad],
       py_modules = ["sparsepad"],
       )