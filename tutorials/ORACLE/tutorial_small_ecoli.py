# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2018 Laboratory of Computational Systems Biotechnology (LCSB),
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
from pytfa.io.viz import get_reaction_data

from skimpy.core import *
from skimpy.mechanisms import *
from skimpy.utils.namespace import *
from skimpy.sampling.simple_parameter_sampler import SimpleParameterSampler
from skimpy.core.solution import ODESolutionPopulation
from skimpy.io.generate_from_pytfa import FromPyTFA
from skimpy.utils.general import sanitize_cobra_vars
from skimpy.utils.tabdict import TabDict

from skimpy.analysis.oracle import *

CPLEX = 'optlang-cplex'
GLPK = 'optlang-glpk'

""" 
Import and curate a model
"""

cobra_model = import_matlab_model('../../models/small_ecoli.mat')
#Test the model
solution = cobra_model.optimize()

""" 
Make tfa analysis of the model
"""
thermo_data = load_thermoDB('../../data/small_ecoli.thermodb')

tmodel= pytfa.ThermoModel(thermo_data, cobra_model)
# for comp in tmodel.compartments.values():
#     comp['c_min'] = 1e-8

tmodel.prepare()
tmodel.convert(add_displacement = True)

# Set the solver
tmodel.solver = CPLEX
# Set solver options
tmodel.solver.configuration.tolerances.optimality = 1e-9
tmodel.solver.configuration.tolerances.feasibility = 1e-9
tmodel.solver.configuration.tolerances.integrality = 1e-9

# Find a solution
solution = tmodel.optimize()


"""
Prepare the model to sample parameters 
"""

# Add minimum flux requirements basal fluxes 1e-6
# safe: ensure that fluxes that cant obey the minimum requirement are removed
basal_flux = 1e-6 # mmol/gDW/hr
tmodel = add_min_flux_requirements(tmodel, basal_flux, inplace=True)
solution = tmodel.optimize()

# Fix the flux directionality profile (FDP)
tmodel = fix_directionality(tmodel, solution, inplace=True)
solution = tmodel.optimize()

# Add dummy free energy constrains for reaction of unknown free energy
tmodel = add_undefined_delta_g(tmodel, solution, delta_g_std=0.0, delta_g_std_err=100000.0, inplace=True)
solution = tmodel.optimize()

# Force a minimal thermodynamic displacement
min_log_displacement = 1e-2
tmodel = add_min_log_displacement(tmodel,min_log_displacement)
solution = tmodel.optimize()


"""
Sample the concentration and flux space for the FDP 
to find reference fluxes and concentrations
"""


"""
Convert fluxes to importable format 
"""

# Map fluxes back to reaction variables
this_flux_solution = get_reaction_data(tmodel, solution.raw)
# Create the flux dict
# Convert fluxes from mmol/gDW/hr to mol/L/s
# eColi 0.39 gDW/L
flux_dict = (0.39*1e-3*this_flux_solution[[i.id for i in tmodel.reactions]]).to_dict()

# Create a concentration dict with consistent names
variable_names = tmodel.log_concentration.list_attr('name')
metabolite_ids = tmodel.log_concentration.list_attr('id')
#Get conentrations in mol
temp_concentration_dict = np.exp(solution.raw[variable_names]).to_dict()

# Map concentration names
mapping_dict = {k:sanitize_cobra_vars(v) for k,v in zip(variable_names,metabolite_ids)}
concentration_dict = {mapping_dict[k]:v for k,v in temp_concentration_dict.items()}

""" 
Import the Kinetic Model 
"""

small_molecules = ['h_c','h_e','h_m',
                   'h2o2_c','h2o2_e',
                   'co2_c','co2_r','co2_e',' co2_m',
                   'pi_c','pi_r','pi_e','pi_m',
                   'o2_c','o2_r','o2_e','o2_m',
                   'o2s_c', 'o2s_m', 'o2s_e',
                   'ppi_c','ppi_m','ppi_r',
                   'hco3_c','hco3_e','hco3_m',
                   'na1_c','na1_e']

model_gen = FromPyTFA(reactants_to_exclude = ['h2o_e','h2o_c'])
kmodel = model_gen.import_model(tmodel, solution.raw)


"""
Compile the model
"""
kmodel.prepare(mca=True)
# Compile MCA functions
kmodel.compile_mca(sim_type=QSSA)


"""
Sample the kinetic parameters using the log linear formulation
of the Jacobian as described in: 
Wang, L. Q., I. Birol and V. Hatzimanikatis (2004). "Metabolic control analysis under uncertainty:
Framework development and case studies." Biophysical Journal 87(6): 3750-3763.
"""

# Initialize parameter sampler
sampling_parameters = SimpleParameterSampler.Parameters(n_samples=100)
sampler = SimpleParameterSampler(sampling_parameters)

# Sample the model
parameter_population = sampler.sample(kmodel, flux_dict, concentration_dict)



"""
Integrate the ODEs
"""

kmodel.compile_ode(sim_type=QSSA)
kmodel.initial_conditions = TabDict([(k,v)for k,v in concentration_dict.items()])

solutions = []
for parameters in parameter_population:
    kmodel.ode_fun.parameter_values = parameters
    this_sol_qssa = kmodel.solve_ode(np.linspace(0.0, 10.0, 1000), solver_type='cvode')
    solutions.append(this_sol_qssa)

this_sol_qssa.plot('output/tutorial_oracle_small_ecoli.html')

solpop = ODESolutionPopulation(solutions)
#If desired plot the solution population for every metabolite
#solpop.plot('output/tutorial_oracle_pop_{}.html')
