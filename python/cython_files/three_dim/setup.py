# python f-setup.py build_ext --inplace
#   cython f.pyx -> f.cpp
#   g++ -c f.cpp -> f.o
#   g++ -c fc.cpp -> fc.o
#   link f.o fc.o -> f.so

# distutils uses the Makefile distutils.sysconfig.get_makefile_filename()
# for compiling and linking: a sea of options.

# http://docs.python.org/distutils/introduction.html
# http://docs.python.org/distutils/apiref.html  20 pages ...
# https://stackoverflow.com/questions/tagged/distutils+python

import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import sys
# from Cython.Build import cythonize

ext_modules = [
    Extension(
        name="wmm3D_cython",
        sources=["wrappers/wmm3D_cython.pyx", "wrappers/wmm3D_c.cpp"],
        compiler_directives={'language_level' : sys.version_info[0]},
        include_dirs = [numpy.get_include()],
        extra_compile_args=['-std=c++14'],
        language="c++",
    ),
    Extension(
        name="aug_wmm3D_cython",
        sources=["wrappers/aug_wmm3D_cython.pyx", "wrappers/aug_wmm3D_c.cpp"],
        compiler_directives={'language_level' : sys.version_info[0]},
        include_dirs = [numpy.get_include()],
        extra_compile_args=['-std=c++14'],
        language="c++",
    ),
    Extension(
        name="ani_wmm3D_cython",
        sources=["wrappers/ani_wmm3D_cython.pyx", "wrappers/ani_wmm3D_c.cpp"],
        compiler_directives={'language_level' : sys.version_info[0]},
        include_dirs = [numpy.get_include()],
        extra_compile_args=['-std=c++14'],
        language="c++",
    ),
    # Extension(
    #     name="ani_wmm3D_riemann_cython",
    #     sources=["wrappers/ani_wmm3D_riemann_cython.pyx", "wrappers/ani_wmm3D_riemann_c.cpp"],
    #     compiler_directives={'language_level' : sys.version_info[0]},
    #     include_dirs = [numpy.get_include()],
    #     extra_compile_args=['-std=c++14'],
    #     language="c++",
    # ),
    Extension(
        name="aug_ani_wmm3D_cython",
        sources=["wrappers/aug_ani_wmm3D_cython.pyx", "wrappers/aug_ani_wmm3D_c.cpp"],
        compiler_directives={'language_level' : sys.version_info[0]},
        include_dirs = [numpy.get_include()],
        extra_compile_args=['-std=c++14'],
        language="c++",
    )
]

setup(
    name = 'mm',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)