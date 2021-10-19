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
from skimpy.mechanisms.generalized_reversible_hill_n_n_h1_with_inhibition import *


# Make the PFK a generalized reversible Hill with inhibition
name = 'pfk'

SpecificHill = make_generalized_reversible_hill_n_n_h1_with_inhibition([-1, -1, 1, 1], [1, 1])
metabolites = SpecificHill.Reactants(substrate1 = 'A',
                                     substrate2 = 'B',
                                     product1 = 'C',
                                     product2 = 'D')

inhibitors = SpecificHill.Inhibitors(inhibitor1 = 'I',
                                     inhibitor2 = 'J')


# QSSA Method
parameters = SpecificHill.Parameters(
                            vmax_forward = 1.0,
                            k_equilibrium=2.0,
                            km_substrate1 = 10.0,
                            km_substrate2 = 10.0,
                            km_product1 = 10.0,
                            km_product2 = 10.0,
                            ki_inhibitor1 = 1.0,
                            ki_inhibitor2 = 1.0)

pfk = Reaction(name=name,
               mechanism=SpecificHill,
               reactants=metabolites,
               inhibitors=inhibitors,
               )



name = 'inhib1'
metabolites = ReversibleMichaelisMenten.Reactants(substrate = 'C',
                                                   product = 'I')

## QSSA Method
parameters_inh1 = ReversibleMichaelisMenten.Parameters(
    vmax_forward = 1.0,
    k_equilibrium=2.0,
    km_substrate = 10.0,
    km_product = 10.0,
    total_enzyme_concentration = 1.0,
)


inhib1 = Reaction(name=name,
               mechanism=ReversibleMichaelisMenten,
               reactants=metabolites,
               )

"""
reactants = make_reactant_set(mechanism,
                                  this_reaction_reactants,
                                  reactant_relations)
"""

name = 'inhib2'
metabolites = ReversibleMichaelisMenten.Reactants(substrate = 'D',
                                                   product = 'J')

## QSSA Method
parameters_inh2 = ReversibleMichaelisMenten.Parameters(
    vmax_forward = 1.0,
    k_equilibrium=2.0,
    km_substrate = 10.0,
    km_product = 10.0,
    total_enzyme_concentration = 1.0,
)


inhib2 = Reaction(name=name,
               mechanism=ReversibleMichaelisMenten,
               reactants=metabolites,
               )

this_model = KineticModel()
this_model.add_reaction(pfk)
this_model.add_reaction(inhib1)
this_model.add_reaction(inhib2)
this_model.parametrize_by_reaction({inhib1.name:parameters_inh1,
                                    inhib2.name:parameters_inh2,
                                    pfk.name: parameters})

this_model.compile_ode(sim_type = QSSA)

this_model.initial_conditions['A'] = 10.0
this_model.initial_conditions['B'] = 7.0
this_model.initial_conditions['C'] = 0.0
this_model.initial_conditions['D'] = 0.25
this_model.initial_conditions['I'] = 0.0
this_model.initial_conditions['J'] = 6.0

this_sol_qssa = this_model.solve_ode(np.linspace(0.0, 50.0, 500),solver_type = 'cvode')

this_sol_qssa.plot('./output/base_out_qssa_hill.html')
