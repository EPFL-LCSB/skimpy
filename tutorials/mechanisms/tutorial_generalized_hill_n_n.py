# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2017 Laboratory of Computational Systems Biotechnology (LCSB),
Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""

import numpy as np
# Test models
from skimpy.core import *
from skimpy.mechanisms import *

name = 'pfk'

SpecificHill= make_generalized_reversible_hill_n_n([-1, -1, 1, 1])
metabolites = SpecificHill.Reactants(substrate1 = 'A',
                                     substrate2 = 'B',
                                     product1 = 'C',
                                     product2 = 'D')

# thermo_data = {'S':     1e-2,
#                'P':     1e-2,
#                'sig_S': 0.1,
#                'sig_P': 0.1,
#                'gamma': 0.1,
#                'flux':  1.0,
#                'E_tot': 1e-5}

## QSSA Method
parameters = SpecificHill.Parameters(
                            vmax_forward = 1.0,
                            k_equilibrium=2.0,
                            km_substrate1 = 10.0,
                            km_substrate2 = 10.0,
                            km_product1 = 10.0,
                            km_product2 = 10.0,
                            hill_coefficient = 1.0)

pfk = Reaction(name=name,
               mechanism=SpecificHill,
               reactants=metabolites,
               )

this_model = KineticModel()
this_model.add_reaction(pfk)
this_model.parametrize_by_reaction({pfk.name:parameters})
this_model.compile_ode(sim_type = QSSA)

this_model.initial_conditions['A'] = 10.0
this_model.initial_conditions['B'] = 10.0
this_model.initial_conditions['C'] = 10.0
this_model.initial_conditions['D'] = 10.0


this_sol_qssa = this_model.solve_ode(np.linspace(0.0, 50.0, 500),solver_type = 'cvode')

this_sol_qssa.plot('output/gen_hill_out_qssa.html')