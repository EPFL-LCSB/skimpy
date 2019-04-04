SKiMpy: Symbolic Kinetic Models in Python
==========================================

This package implements a toolbox to build and analyse Kinetic Models.
The toolbox uses symbolic expressions to build the ODEs that represent the kinetic model. Within the package the following
methods are implemented:

    - Steady state consistent parameter sampling(see ORACLE)
    - Metabolic control analysis
    - Modal analysis

Requirements
------------

This module was developed in Python 3.5, and it is recommended to run Python 3.5.
The module also was tested in Python 3.6.


Further the following pip-python packages are required
    - sympy >= 1.1.
    - pytest
    - scipy
    - numpy
    - bokeh
    - pandas
    - Cython
    - scikits.odes
    - deap

