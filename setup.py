#!/usr/bin/env python3

""" Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team


"""

from setuptools import setup, find_packages
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

from numpy.distutils.system_info import default_include_dirs, default_lib_dirs

from distutils.sysconfig import get_config_vars

import os
import sys


if sys.platform == 'win32':
    libraries = ["gmp", "flint", "mpir", "mpfr", "pthreads"]
else:
    libraries = ["gmp", "flint"]
    (opt,) = get_config_vars('OPT')
    os.environ['OPT'] = " ".join(flag for flag in opt.split() if flag != '-Wstrict-prototypes')

default_include_dirs += [
    os.path.join(d, "flint") for d in default_include_dirs
]

# Current version 0.9,1-beta
version_tag = '0.9.1-beta'


extensions = [
            Extension("skimpy.nullspace",
                      sources=["skimpy/cython/nullspace.pyx"],
                      language = 'c',
                      libraries=libraries,
                      library_dirs = default_lib_dirs,
                      include_dirs = default_include_dirs),
            ]

setup(name='skimpy',
      version=version_tag,
      author='SKiMPy team',
      author_email='softwares.lcsb@epfl.ch',
      url='https://github.com/EPFL-LCSB/skimpy/',
      download_url='https://github.com/EPFL-LCSB/skimpy/archive/'+version_tag+'.tar.gz',
      install_requires=['sympy',
                        'pytest',
                        'scipy',
                        'numpy',
                        'pandas',
                        'Cython',
                        'markupsafe<=2.0.1',
                        'bokeh>=0.12.0',
                        'scikits.odes==2.6.3',
                        'deap',
                        'dill',
                        'h5py',
                        'escher',
                        'matplotlib',
                        'pytfa',
                        ],
      packages = find_packages(),
      python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4',
      description='SKiMPy allows to develop and analyze large-scale biochemical kinetic models',
      keywords=['skimpy', 'kinetic', 'models'],

      #ext_modules=cythonize(extensions ),
      ext_modules=extensions,
      cmdclass={'build_ext': build_ext},

      license='Apache2',

      # See https://PyPI.python.org/PyPI?%3Aaction=list_classifiers
      classifiers=[
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            'Development Status :: 3 - Alpha',

            # Indicate who your project is intended for
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry'
            'Environment :: Console',

            # Pick your license as you wish (should match "license" above)
            'License :: OSI Approved :: Apache Software License',

            # Specify the Python versions you support here. In particular, ensure
            # that you indicate whether you support Python 2, Python 3 or both.
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',

      ],
     )

