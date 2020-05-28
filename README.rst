SKiMpy: Symbolic Kinetic Models in Python
==========================================
|Build Status| |Codecov| |Codacy branch grade| |license| 

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

.. |license| image:: http://img.shields.io/badge/license-APACHE2-blue.svg
   :target: https://github.com/EPFL-LCSB/skimpy/blob/master/LICENSE.txt
.. |Build Status| image:: https://travis-ci.org/EPFL-LCSB/skimpy.svg?branch=master
   :target: https://travis-ci.org/EPFL-LCSB/skimpy
.. |Codecov| image:: https://img.shields.io/codecov/c/github/EPFL-LCSB/skimpy.svg
   :target: https://codecov.io/gh/EPFL-LCSB/skimpy
.. |Codacy branch grade| image:: https://img.shields.io/codacy/grade/d56d598a8a3b444e8ea5fb1f7eee6e2a
   :target: https://www.codacy.com/app/realLCSB/skimpy
