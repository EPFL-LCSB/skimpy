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
import numpy as np
from skimpy.core import *
from skimpy.mechanisms import *

name = 'UniBi'
metabolites = UniBiReversibleHill.Reactants(substrate='A',
                                            product1='B',
                                            product2='C',
                                            )

## QSSA Method
parameters = UniBiReversibleHill.Parameters(
    vmax_forward=1.0,
    k_equilibrium=2.0,
    hill_coefficient=1,
    km_substrate=10.0,
    km_product1=10.0,
    km_product2=10.0,
)


pfk = Reaction(name=name,
               mechanism = UniBiReversibleHill,
               reactants=metabolites,
               )

this_model = KineticModel()
this_model.add_reaction(pfk)
this_model.parametrize_by_reaction({pfk.name:parameters})
this_model.compile_ode(sim_type = QSSA)

this_model.initial_conditions['A'] = 1.0
this_model.initial_conditions['B'] = 1.0

this_sol_qssa = this_model.solve_ode(np.linspace(0.0, 100.0, 1000), solver_type='cvode')

this_sol_qssa.plot('output/uni_bi_base_out_qssa.html')




"""
BiBi Michaelis Menten Kinetics
"""
name = 'BiUni'
metabolites = BiUniReversibleHill.Reactants(substrate1 = 'A',
                                            substrate2 = 'B',
                                            product = 'C',
                                            )

parameters = BiUniReversibleHill.Parameters(
    vmax_forward=1.0,
    k_equilibrium=5.0,
    hill_coefficient=1,
    km_substrate1=10,
    km_substrate2=10,
    km_product=1.0,
)

hxk = Reaction(name=name,
               mechanism=BiUniReversibleHill,
               reactants=metabolites,
               )

this_model = KineticModel()
this_model.add_reaction(hxk)
this_model.parametrize_by_reaction({hxk.name:parameters})
this_model.compile_ode(sim_type = QSSA)

this_model.initial_conditions['A'] = 100.0
this_model.initial_conditions['B'] = 1.0
this_model.initial_conditions['C1'] = 3.0
this_model.initial_conditions['C2'] = 5.0

this_sol_qssa = this_model.solve_ode(np.linspace(0.0, 10.0, 1000),solver_type = 'cvode')

this_sol_qssa.plot('output/bi_uni_base_out_qssa.html')
