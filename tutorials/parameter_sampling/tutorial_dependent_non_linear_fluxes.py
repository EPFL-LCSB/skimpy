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
from skimpy.sampling.simple_parameter_sampler import SimpleParameterSampler
from skimpy.core.solution import ODESolutionPopulation

# 'DM_13dpg'        10
# 'DM_2h3oppan'     -2
# 'DM_pep'          -8
# 'DM_atp'         -12
# 'DM_adp'          12
# 'DM_nadh'         -2
# 'DM_nad'           2
# 'DM_h'           -12
# 'DM_h2o'          -8
# 'PGK'             10
# 'PGK'             10
# 'GLYCK'            9
# 'PGM'              1
# 'GLYCK2'           7
# 'TRSARr'           2
# 'ENO'              8

flux_dict = {'PGK'         : 10.0,
             'PGM'         : 1.0,
             'GLYCK'       : 9.0,
             'GLYCK2'      : 7.0,
             'TRSARr'      : 2.0,
             'ENO'         : 8.0,
             'DM_13dpg'    : 10,
             'DM_2h3oppan' : -2,
             'DM_pep'      : -8,
             'DM_atp'      : -12,
             'DM_adp'      : 12,
             'DM_nadh'     : -2,
             'DM_nad'      : 2,
             }

# Build linear Pathway model
metabolites_pgk    = RandBiBiReversibleMichaelisMenten.Reactants(
    substrate1     = '_13dpg',
    substrate2     = 'adp',
    product1       = '_3pg',
    product2       = 'atp')
metabolites_pgm    = ReversibleMichaelisMenten.Reactants(
    substrate      = '_3pg',
    product        = '_2pg')
metabolites_glyck  = RandBiBiReversibleMichaelisMenten.Reactants(
    substrate1     = '_3pg',
    substrate2     = 'adp',
    product1       = 'glyc',
    product2       = 'atp')
metabolites_glyck2 = RandBiBiReversibleMichaelisMenten.Reactants(
    substrate1     = 'glyc',
    substrate2     = 'atp',
    product1       = '_2pg',
    product2       = 'adp')
metabolites_trsarr = RandBiBiReversibleMichaelisMenten.Reactants(
    substrate1     = 'glyc',
    substrate2     = 'nad',
    product1       = '_2h3oppan',
    product2       = 'nadh')
metabolites_eno    = ReversibleMichaelisMenten.Reactants(
    substrate      = '_2pg',
    product        = 'pep')


def keq(deltag):
    return np.exp(-deltag*4184/(8.31*298.15))

## QSSA Method
parameters_pgm    = ReversibleMichaelisMenten.Parameters(k_equilibrium=keq(-1.0))
parameters_pgk    = RandBiBiReversibleMichaelisMenten.Parameters(k_equilibrium=keq(-0.7996))
parameters_glyck  = RandBiBiReversibleMichaelisMenten.Parameters(k_equilibrium=keq(-5.4833))
parameters_glyck2 = RandBiBiReversibleMichaelisMenten.Parameters(k_equilibrium=keq(-5.3603))
parameters_trsarr = RandBiBiReversibleMichaelisMenten.Parameters(k_equilibrium=keq(-3.7509))
parameters_eno    = ReversibleMichaelisMenten.Parameters(k_equilibrium=keq(-1.6601))

pgk    = Reaction(name='PGK',
                  mechanism=RandBiBiReversibleMichaelisMenten,
                  reactants=metabolites_pgk,
                 )
pgm    = Reaction(name='PGM',
                  mechanism=ReversibleMichaelisMenten,
                  reactants=metabolites_pgm,
                  )
glyck  = Reaction(name='GLYCK',
                  mechanism=RandBiBiReversibleMichaelisMenten,
                  reactants=metabolites_glyck,
                  )
glyck2 = Reaction(name='GLYCK2',
                  mechanism=RandBiBiReversibleMichaelisMenten,
                  reactants=metabolites_glyck2,
                  )
trsarr = Reaction(name='TRSARr',
                  mechanism=RandBiBiReversibleMichaelisMenten,
                  reactants=metabolites_trsarr,
                  )
eno    = Reaction(name='ENO',
                  mechanism=ReversibleMichaelisMenten,
                  reactants=metabolites_eno,
                  )


this_model = KineticModel()
this_model.add_reaction(pgk)
this_model.add_reaction(pgm)
this_model.add_reaction(glyck)
this_model.add_reaction(glyck2)
this_model.add_reaction(trsarr)
this_model.add_reaction(eno)

the_boundary_condition = BoundaryFlux("_13dpg", flux_dict['DM_13dpg'])
this_model.add_boundary_condition(the_boundary_condition)

the_boundary_condition = BoundaryFlux("pep", flux_dict['DM_pep'])
this_model.add_boundary_condition(the_boundary_condition)

the_boundary_condition = BoundaryFlux("_2h3oppan", flux_dict['DM_2h3oppan'])
this_model.add_boundary_condition(the_boundary_condition)

the_boundary_condition = BoundaryFlux("atp", flux_dict['DM_atp'])
this_model.add_boundary_condition(the_boundary_condition)

the_boundary_condition = BoundaryFlux("adp", flux_dict['DM_adp'])
this_model.add_boundary_condition(the_boundary_condition)

the_boundary_condition = BoundaryFlux("nad", flux_dict['DM_nad'])
this_model.add_boundary_condition(the_boundary_condition)

the_boundary_condition = BoundaryFlux("nadh", flux_dict['DM_nadh'])
this_model.add_boundary_condition(the_boundary_condition)


this_model.parametrize_by_reaction({'PGK'   : parameters_pgk,
                        'PGM'   : parameters_pgm,
                        'GLYCK' : parameters_glyck,
                        'GLYCK2': parameters_glyck2,
                        'TRSARr': parameters_trsarr,
                        'ENO'   : parameters_eno})


this_model.prepare(mca=True)
this_model.compile_mca()


concentration_dict = {'_13dpg'    : 1.0,
                      '_2pg'      : 1.0,
                      '_3pg'      : 1.0,
                      'glyc'      : 1.0,
                      'pep'       : 1.0,
                      '_2h3oppan' : 1.0,
                      'atp'       : 1.0,
                      'adp'       : 1.0,
                      'nad'       : 1.0,
                      'nadh'      : 1.0,
                      }


parameters = SimpleParameterSampler.Parameters(n_samples=10)
sampler = SimpleParameterSampler(parameters)

parameter_population = sampler.sample(this_model, flux_dict, concentration_dict)

this_model.compile_ode(sim_type = QSSA)

#
this_model.initial_conditions['_13dpg']     = 5.0
this_model.initial_conditions['_2pg']       = 1.0
this_model.initial_conditions['_3pg']       = 1.0
this_model.initial_conditions['glyc']       = 1.0
this_model.initial_conditions['pep']        = 1.0
this_model.initial_conditions['_2h3oppan']  = 1.0
this_model.initial_conditions['atp']        = 1.0
this_model.initial_conditions['adp']        = 2.0
this_model.initial_conditions['nad']        = 1.0
this_model.initial_conditions['nadh']       = 1.0

this_model.logger.setLevel('INFO')

solutions = []
for parameters in parameter_population:
    this_model.ode_fun.parameter_values = parameters
    this_sol_qssa = this_model.solve_ode(np.linspace(0.0, 50.0, 500), solver_type='cvode')
    solutions.append(this_sol_qssa)

this_sol_qssa.plot('output/non_linear_qssa.html')


solpop = ODESolutionPopulation(solutions)
solpop.plot('output/non_linear_ode_pop_{}.html')

