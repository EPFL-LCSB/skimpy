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
from skimpy.core.solution import ODESolutionPopulation
from skimpy.sampling.simple_parameter_sampler import SimpleParameterSampler

# Build linear Pathway model
metabolites_1 = ReversibleMichaelisMenten.Reactants(substrate='A', product='B')
metabolites_2 = ReversibleMichaelisMenten.Reactants(substrate='B', product='C')
metabolites_3 = ReversibleMichaelisMenten.Reactants(substrate='C', product='D')

## QSSA Method
parameters_1 = ReversibleMichaelisMenten.Parameters(k_equilibrium=1.5)
parameters_2 = ReversibleMichaelisMenten.Parameters(k_equilibrium=2.0)
parameters_3 = ReversibleMichaelisMenten.Parameters(k_equilibrium=3.0)

reaction1 = Reaction(name='E1',
                     mechanism=ReversibleMichaelisMenten,
                     reactants=metabolites_1,
                     )

reaction2 = Reaction(name='E2',
                     mechanism=ReversibleMichaelisMenten,
                     reactants=metabolites_2,
                     )

reaction3 = Reaction(name='E3',
                     mechanism=ReversibleMichaelisMenten,
                     reactants=metabolites_3,
                     )

this_model = KineticModel()
this_model.add_reaction(reaction1)
this_model.add_reaction(reaction2)
this_model.add_reaction(reaction3)

the_boundary_condition = ConstantConcentration(this_model.reactants['A'])
this_model.add_boundary_condition(the_boundary_condition)

the_boundary_condition = ConstantConcentration(this_model.reactants['D'])
this_model.add_boundary_condition(the_boundary_condition)

this_model.parametrize_by_reaction({'E1': parameters_1,
                                    'E2': parameters_2,
                                    'E3': parameters_3})

this_model.prepare(mca=True)
#this_model.compile_mca()

this_model.compile_jacobian(type=SYMBOLIC,sim_type = QSSA)

flux_dict = {'E1': 1.0, 'E2': 1.0, 'E3': 1.0}
concentration_dict = {'A': 3.0, 'B': 2.0, 'C': 1.0, 'D': 0.5}


parameters = SimpleParameterSampler.Parameters(n_samples=10)
sampler = SimpleParameterSampler(parameters)

parameter_population = sampler.sample(this_model, flux_dict,
                                      concentration_dict)

#this_model.compile_ode(sim_type = QSSA)

# TODO can change this back to this_model.initial_conditions.A = 3.0 once
# tabdict is fixed
this_model.initial_conditions['A'] = 3.0
this_model.initial_conditions['B'] = 2.0
this_model.initial_conditions['C'] = 3.0
this_model.initial_conditions['D'] = 0.5

solutions = []
for parameters in parameter_population:
    # TODO system is too small for sparse eig handle properly
    this_model.parameters = parameters
    this_sol_qssa = this_model.solve_ode(np.linspace(0.0, 50.0, 500),
                                         solver_type='cvode')
    solutions.append(this_sol_qssa)

this_sol_qssa.plot('output/non_linear_qssa.html')

solpop = ODESolutionPopulation(solutions)
solpop.plot('output/non_linea_pop_{}.html')
