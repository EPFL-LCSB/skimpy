SKiMpy: Symbolic Kinetic Models in Python
==========================================
 |docs| |Build Status| |Codecov| |Codacy branch grade| |license| 

SKiMpy is a python package bridging implementing an efficient kinetic model-ing toolbox to build and analyze large-scale kinetic models for various biological domains such as signaling, gene expression, and metabolism. Furthermore, we demonstrate how this toolbox is used to parameterize kinetic models around a steady-state reference efficiently. Finally, we show how SKiMpy can implement multispecies bioreactor simulations to assess biotechnological processes.


    - Non-linear ordinary equations for large scale kinetic models
    - Steady state consistent parameter sampling(see ORACLE)
    - Metabolic control analysis
    - Modal analysis
    - Uncertainty propagation in metabolic control
    - Multispecies bioreactor modeling

Requirements
------------

You will need to have `Git-LFS <https://git-lfs.github.com/>`_ in order to properly download some binary files:

.. code:: bash

    git clone https://github.com/EPFL-LCSB/skimpy.git /path/to/skimpy
    cd /path/to/skimpy
    git lfs install
    git lfs pull
 


This module was developed in Python 3.6, and it is recommended to run Python 3.6.
The module also was tested in Python 3.7.

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
    - dill
    - h5py
    - escher

Container-based install
-----------------------

You might want to use this program inside of a container. The
|docker|_
subfolder has all the necessary information and source files to set it
up.

.. |docker| replace:: ``docker/``
.. _docker: https://github.com/EPFL-LCSB/skimpy/tree/master/docker


Setup
=====

*This step is not required if you're using the container, which bundles all this.*

If you prefer a local installation of skimpy we provide a conda package bundeling binary depdendecies for linux and OSX operatings systems. 
Windows users can install the conda packages using a linux subsystem via `WLS <https://docs.microsoft.com/en-us/windows/wsl/install>`_.

A comprehensive explanation to install anaconda/miniconda within the WLS linux subsystem can be found `here <https://gist.github.com/kauffmanes/5e74916617f9993bc3479f401dfec7da>`_.

.. code:: bash

  conda TO BE FINISHED SOON 



Alternatively you can install this module from source using ``pip``:
*For Python 3, you might have to use* ``pip3`` *instead of* ``pip``

.. code:: bash

    git clone https://github.com/EPFL-LCSB/skimpy.git /path/to/skimpy
    pip3 install -e /path/to/skimpy
    
    
Note that the manual installation requieres additional packages:
  - cython3
  - build-essential
  - gcc
  - gfortran
  - libsundials-dev
  - libflint- dev


Installation of these packages on the native windows system is challenging we thus recomend windows users to install the package only with in linux subsystem using `WLS <https://docs.microsoft.com/en-us/windows/wsl/install>`_.
  
To use the ODE integration features scikit-odes is required to be installed beforehand follwoing the instructions found `here <https://scikits-odes.readthedocs.io/en/stable/installation.html>`_. To use the 'cvode' solver from the scikit-odes packages we strongly recomend to install the `sundials <https://computing.llnl.gov/projects/sundials>`_ solvers as ODE integration of large ODE system can be quite slow with python implemented solvers see benchmark `here <https://scikits-odes.readthedocs.io/en/stable/solvers.html>`_.

Installation from source has been tested on Ubuntu 21.10  (`@eladnoor <https://github.com/eladnoor/>`_) the additional packages can be installed using:

.. code:: bash

  sudo apt install cython3 gfortran libsundials-dev libflint-dev
 
Windows users using *WSL* can install these dependicies in a similar fashion after starting the subsystem console. 


Quick start
===========
To get right into building kinetic models please find below a simple example to get started:

.. code-block:: python

    import numpy as np
    from skimpy.core import *
    from skimpy.mechanisms import *

    name = 'pfk'
    metabolites = ReversibleMichaelisMenten.Reactants(substrate = 'A',
                                                       product = 'B')

    parameters = ReversibleMichaelisMenten.Parameters(
        vmax_forward = 1.0,
        k_equilibrium=2.0,
        km_substrate = 10.0,
        km_product = 10.0,
        total_enzyme_concentration = 1.0,
    )


    pfk = Reaction(name=name,
                   mechanism = ReversibleMichaelisMenten,
                   reactants=metabolites,
                   )

    this_model = KineticModel()
    this_model.add_reaction(pfk)
    this_model.parametrize_by_reaction({pfk.name:parameters})
    this_model.compile_ode(sim_type = QSSA)

    this_model.initial_conditions['A'] = 1.0
    this_model.initial_conditions['B'] = 1.0

    this_sol_qssa = this_model.solve_ode(np.linspace(0.0, 100.0, 1000), solver_type='cvode')

    this_sol_qssa.plot('output/uni_uni_base_out_qssa.html')


More information can be found
`here <http://real-skimpy.readthedocs.io/en/latest/quickstart.html>`__.


   
License
========

The software in this repository is put under an APACHE-2.0 licensing scheme - please see the `LICENSE <https://github.com/EPFL-LCSB/skimpy/blob/master/LICENSE.txt>`_ file for more details.

.. |docs| image:: https://readthedocs.org/projects/real-skimpy/badge/?version=latest
   :target: https://real-skimpy.readthedocs.io/en/latest/?badge=latest
.. |license| image:: http://img.shields.io/badge/license-APACHE2-blue.svg
   :target: https://github.com/EPFL-LCSB/skimpy/blob/master/LICENSE.txt
.. |Build Status| image:: https://travis-ci.org/EPFL-LCSB/skimpy.svg?branch=master
   :target: https://travis-ci.org/EPFL-LCSB/skimpy
.. |Codecov| image:: https://img.shields.io/codecov/c/github/EPFL-LCSB/skimpy.svg
   :target: https://codecov.io/gh/EPFL-LCSB/skimpy
.. |Codacy branch grade| image:: https://img.shields.io/codacy/grade/d56d598a8a3b444e8ea5fb1f7eee6e2a
   :target: https://www.codacy.com/app/realLCSB/skimpy
