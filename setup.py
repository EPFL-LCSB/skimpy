""" Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team


"""

from setuptools import setup, find_packages
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext




version_tag = '0.0.1'


extensions = [
            Extension("skimpy.nullspace",
                      sources=["skimpy/cython/nullspace.pyx"],
                      library_dirs=['/usr/lib'],
                      extra_compile_args=[],
                      extra_link_args=['-lgmp','-lflint'],
                      include_dirs=[],
                      language = 'c',
                      libraries=['/usr/local/lib/libgmp.so.23',
                                 '/usr/local/lib/libflint.so.13'],
                      )
            ,]


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
                        'bokeh',
                        'Cython',
                        'scikits.odes==2.4.1',
                        'deap',
                        'dill',
                        'h5py',
                        'escher',
                        'matplotlib',
                        ],
      packages = find_packages(),
      python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4',
      description='SKiMPy adds Thermodynamics-based Flux Analysis',
      keywords=['skimpy','kinetic','models'],
      extras_require={ 'ORACLE':  ["pytfa"], },

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
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
      ],
     )
