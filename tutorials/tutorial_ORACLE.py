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

import pytfa
from pytfa.io import import_matlab_model, load_thermoDB
from pytfa.optim.variables import DeltaG,DeltaGstd,ThermoDisplacement
from pytfa.analysis import  variability_analysis,           \
                            apply_reaction_variability,     \
                            apply_generic_variability,       \
                            apply_directionality

from skimpy.core import *
from skimpy.mechanisms import *
from skimpy.utils.namespace import *

from skimpy.sampling import SimpleParameterSampler
from skimpy.core.solution import ODESolutionPopulation
from skimpy.io.generate_from_pytfa import FromPyTFA
from skimpy.utils.general import sanitize_cobra_vars
from skimpy.utils.tabdict import TabDict

""" Import an curate the model """
#this_cobra_model = import_matlab_model('../models/toy_model.mat', 'ToyModel_DP')
this_cobra_model = import_matlab_model('../models/toy_model_maria.mat', 'model')


""" Make tfa analysis of the model """

# Convert to a thermodynamics model
thermo_data = load_thermoDB('../data/thermo_data.thermodb')
this_pytfa_model = pytfa.ThermoModel(thermo_data, this_cobra_model)

GLPK = 'optlang-glpk'
this_pytfa_model.solver = GLPK

## TFA conversion
this_pytfa_model.prepare()
this_pytfa_model.convert(add_displacement=True)
solution = this_pytfa_model.optimize()

min_log_displacement = 1e-3
for ln_gamma in this_pytfa_model.thermo_displacement:
     if ln_gamma.variable.primal > 0:
         ln_gamma.variable.lb = min_log_displacement
     if ln_gamma.variable.primal < 0:
         ln_gamma.variable.ub = -min_log_displacement

solution = this_pytfa_model.optimize()


""" Get a Kinetic Model """
# Generate the KineticModel
# TODO This is really bad we need a better way
small_molecules = ['h_c','h_e']

model_gen = FromPyTFA(water='h2o')
this_skimpy_model = model_gen.  import_model(this_pytfa_model,solution)

# Compile MCA functions
this_skimpy_model.compile_mca(sim_type=QSSA)

# Initialize parameter sampler
sampling_parameters = SimpleParameterSampler.Parameters(n_samples=100)
sampler = SimpleParameterSampler(sampling_parameters)

# Create the flux dict
flux_dict = solution[[i for i in this_skimpy_model.reactions.keys()]].to_dict()

# Create a concentration dict with consistent names
variable_names = this_pytfa_model.log_concentration.list_attr('name')
metabolite_ids = this_pytfa_model.log_concentration.list_attr('id')

temp_concentration_dict = np.exp(solution[variable_names]).to_dict()

# Map concentration names
mapping_dict = {k:sanitize_cobra_vars(v)
                for k,v in zip(variable_names,metabolite_ids)}
concentration_dict = {mapping_dict[k]:v for k,v in temp_concentration_dict.items()}


parameter_population = sampler.sample(this_skimpy_model, flux_dict, concentration_dict)


this_skimpy_model.compile_ode(sim_type=QSSA)
this_skimpy_model.initial_conditions = TabDict([(k,v)for k,v in concentration_dict.items()])

solutions = []
for parameters in parameter_population:
    this_skimpy_model.ode_fun.parameter_values = parameters
    this_sol_qssa = this_skimpy_model.solve_ode(np.linspace(0.0, 10.0, 1000), solver_type='cvode')
    solutions.append(this_sol_qssa)

this_sol_qssa.plot('output/tutorial_oracle.html')

solpop = ODESolutionPopulation(solutions)
solpop.plot('output/tutorial_oracle_pop_{}.html')



