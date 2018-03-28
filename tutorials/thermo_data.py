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
from skimpy.utils.constructors import reversible_michaelis_menten_from_thermo_data

name = 'pfk'
metabolites = ReversibleMichaelisMenten.Substrates(substrate = 'A',
                                                   product = 'B')

thermo_data_pfk = {'concetration_substrate':     1e-2,
                   'concetration_product':       1e-2,
                   'saturation_substrate':       0.1,
                   'saturation_product':         0.1,
                   'thermo_displacement':        0.8,
                   'flux':                       1e-4,
                   'total_enzyme_concentration': 1e-5}

## QSSA Method
parameters = reversible_michaelis_menten_from_thermo_data(thermo_data_pfk)

pfk = Reaction(name=name,
               mechanism = ReversibleMichaelisMenten,
               reactants=metabolites,
               )

this_model = KineticModel()
this_model.add_reaction(pfk)
this_model.parametrize({pfk.name:parameters})
this_model.compile_ode(sim_type = 'QSSA')

this_model.initial_conditions.A = 1e-2
this_model.initial_conditions.B = 1e-2

this_sol_qssa = this_model.solve_ode([0,100.0],solver_type = 'vode')

this_sol_qssa.plot('output/thermo_data_out_qssa.html')

## Full rate method
this_model.compile_ode(sim_type = 'full')

this_model.initial_conditions.A = 1e-2
this_model.initial_conditions.B = 1e-2
this_model.initial_conditions.pfk = (0.8)*thermo_data_pfk['total_enzyme_concentration']
this_model.initial_conditions.EC_pfk = (0.2)*thermo_data_pfk['total_enzyme_concentration']

this_sol_full = this_model.solve_ode([0,100.0], solver_type = 'vode')

this_sol_full.plot('output/thermo_data_out_full.html')


#this_sol_full.species[-1,[2,3]] - this_sol_qssa.species[-1,:]
