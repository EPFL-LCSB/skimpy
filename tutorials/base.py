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

# Test models
from skimpy.core import *
from skimpy.mechanisms import *

name = 'pfk'
metabolites = ReversibleMichaelisMenten.Substrates(substrate = 'A',
                                                   product = 'B')

# thermo_data = {'S':     1e-2,
#                'P':     1e-2,
#                'sig_S': 0.1,
#                'sig_P': 0.1,
#                'gamma': 0.1,
#                'flux':  1.0,
#                'E_tot': 1e-5}

## QSSA Method
parameters = ReversibleMichaelisMenten.Parameters(
    vmax_forward = 1,
    vmax_backward = 0.5,
    km_substrate = 10,
    km_product = 10,
    total_enzyme_concentration = 1,
)

pfk = Reaction(name=name,
               mechanism = ReversibleMichaelisMenten,
               substrates=metabolites,
               )

this_model = KineticModel()
this_model.add_reaction(pfk)
this_model.parametrize({pfk.name:parameters})
this_model.compile_ode(sim_type = 'QSSA')

this_model.initial_conditions.A = 10.0
this_model.initial_conditions.B = 1.0

this_sol_qssa = this_model.solve_ode([0,100.0],solver_type = 'vode')

this_sol_qssa.plot('out_qssa.html')

## Full rate method


this_model.compile_ode(sim_type = 'full')

this_model.initial_conditions.A = 10.0
this_model.initial_conditions.B = 1.0
this_model.initial_conditions.pfk = 1.0

this_sol_full = this_model.solve_ode([0,100.0], solver_type = 'vode')

this_sol_full.plot('out_full.html')
