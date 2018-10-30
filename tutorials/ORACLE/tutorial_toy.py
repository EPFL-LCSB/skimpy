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

from pytfa.optim.variables import ReactionVariable, MetaboliteVariable
from pytfa.io.viz import get_reaction_data

from skimpy.core import *
from skimpy.mechanisms import *
from skimpy.utils.namespace import *
from skimpy.sampling import SimpleParameterSampler
from skimpy.core.solution import ODESolutionPopulation
from skimpy.io.generate_from_pytfa import FromPyTFA
from skimpy.utils.general import sanitize_cobra_vars
from skimpy.utils.tabdict import TabDict

from skimpy.analysis.oracle import *

""" 
Import and curate a model
"""

#this_cobra_model = import_matlab_model('../../models/toy_model.mat', 'ToyModel_DP')
this_cobra_model = import_matlab_model('../../models/toy_model_maria.mat', 'model')


""" 
Make tfa analysis of the model
"""

# Convert to a thermodynamics model
thermo_data = load_thermoDB('../../data/thermo_data.thermodb')
this_pytfa_model = pytfa.ThermoModel(thermo_data, this_cobra_model)

GLPK = 'optlang-glpk'
this_pytfa_model.solver = GLPK

## TFA conversion
this_pytfa_model.prepare()
this_pytfa_model.convert(add_displacement=True)


# We choose a flux directionality profile (FDP)
# with minium fluxes of 1e-3

this_bounds = {  'DM_13dpg': (-10.0, -2.0),
                 'DM_2h3oppan': (1e-3, 100.0),
                 'DM_adp':      (-100.0, -1e-3),
                 'DM_atp':      (1e-3, 100.0),
                 'DM_h':        (1e-3, 100.0),
                 'DM_h2o':      (1e-3, 100.0),
                 'DM_nad':      (-100.0, -1e-3),
                 'DM_nadh':     (1e-3, 100.0),
                 'DM_pep':      (1e-3, 100.0),
                 'ENO':         (2.0, 100.0),
                 'GLYCK':       (1e-3, 100.0),
                 'GLYCK2':      (-100, -1e-3),
                 'PGK':         (1e-3, 100.0),
                 'PGM':         (1e-3, 100.0),
                 'TRSARr':      (2.0, 10.0),
                 'Trp_adp':     (-100.0, 100.0),
                 'Trp_atp':     (-100.0, 100.0),
                 'Trp_h':       (-100.0, 100.0),
                 'Trp_h2o':     (-100.0, 100.0),
                 'Trp_nad':     (-100.0, 100.0),
                 'Trp_nadh':    (-100.0, 100.0)}
# Find a solution for this FDP
solution = this_pytfa_model.optimize()

# Force a minimal thermodynamic displacement
min_log_displacement = 1e-1
add_min_log_displacement(this_pytfa_model,min_log_displacement)

# Find a solution for the model
solution = this_pytfa_model.optimize()


""" 
Get a Kinetic Model 
"""
# Generate the KineticModel

# Define the molecules that should be considered small-molecules
# These molecules will not be accounted explicitly in the kinetic mechanism as
# substrates and products
small_molecules = ['h_c','h_e']

model_gen = FromPyTFA(water='h2o', small_molecules=small_molecules)
this_skimpy_model = model_gen.import_model(this_pytfa_model,solution)

"""
Sanitize the solution to match with the skimpy model
"""

# Map fluxes back to reaction variables
this_flux_solution = get_reaction_data(this_pytfa_model, solution)
# Create the flux dict
flux_dict = this_flux_solution[[i for i in this_skimpy_model.reactions.keys()]].to_dict()

# Create a concentration dict with consistent names
variable_names = this_pytfa_model.log_concentration.list_attr('name')
metabolite_ids = this_pytfa_model.log_concentration.list_attr('id')

temp_concentration_dict = np.exp(solution[variable_names]).to_dict()

# Map concentration names
mapping_dict = {k:sanitize_cobra_vars(v) for k,v in zip(variable_names,metabolite_ids)}
concentration_dict = {mapping_dict[k]:v for k,v in temp_concentration_dict.items()}


"""
Sample the kinetic parameters using MCA
"""

# Compile MCA functions
this_skimpy_model.compile_mca(sim_type=QSSA)

# Initialize parameter sampler
sampling_parameters = SimpleParameterSampler.Parameters(n_samples=100)
sampler = SimpleParameterSampler(sampling_parameters)

parameter_population = sampler.sample(this_skimpy_model, flux_dict, concentration_dict)


"""
Integrate the ODEs
"""

this_skimpy_model.compile_ode(sim_type=QSSA)
this_skimpy_model.initial_conditions = TabDict([(k,v)for k,v in concentration_dict.items()])

solutions = []
for parameters in parameter_population:
    this_skimpy_model.ode_fun.parameter_values = parameters
    this_sol_qssa = this_skimpy_model.solve_ode(np.linspace(0.0, 0.5, 1000), solver_type='cvode')
    solutions.append(this_sol_qssa)

this_sol_qssa.plot('output/tutorial_oracle.html')

solpop = ODESolutionPopulation(solutions)
solpop.plot('output/tutorial_oracle_pop_{}.html')



