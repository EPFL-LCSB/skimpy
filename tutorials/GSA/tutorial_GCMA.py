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
from pytfa.io.viz import get_reaction_data

from skimpy.core import *
from skimpy.mechanisms import *
from skimpy.utils.namespace import *
from skimpy.sampling.simple_resampler import SimpleResampler
from skimpy.sampling.simple_parameter_sampler import SimpleParameterSampler
from skimpy.core.parameters import ParameterValuePopulation
from skimpy.io.generate_from_pytfa import FromPyTFA
from skimpy.utils.general import sanitize_cobra_vars
from skimpy.utils.tabdict import TabDict

from skimpy.analysis.oracle import *

"""
Import and curate a model
"""

this_cobra_model = import_matlab_model('../../models/toy_model.mat',
                                       'model')


"""
Make tfa analysis of the model
"""

# Convert to a thermodynamics model
thermo_data = load_thermoDB('../../data/thermo_data.thermodb')
tmodel = pytfa.ThermoModel(thermo_data, this_cobra_model)

GLPK = 'optlang-glpk'
tmodel.solver = GLPK

# TFA conversion
tmodel.prepare()
tmodel.convert(add_displacement=True)


# We choose a flux directionality profile (FDP)
# with minium fluxes of 1e-3

this_bounds = {'DM_13dpg':    (-10.0, -2.0),
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
solution = tmodel.optimize()

# Force a minimal thermodynamic displacement
min_log_displacement = 1e-1
add_min_log_displacement(tmodel, min_log_displacement)

# Find a solution for the model
solution = tmodel.optimize()


"""
Get a Kinetic Model
"""
# Generate the KineticModel

# Define the molecules that should be considered small-molecules
# These molecules will not be accounted explicitly in the kinetic mechanism as
# substrates and products
small_molecules = ['h_c', 'h_e']

model_gen = FromPyTFA(small_molecules=small_molecules)
kmodel = model_gen.import_model(tmodel, solution.raw)


"""
Sanitize the solution to match with the skimpy model
"""

# Map fluxes back to reaction variables
this_flux_solution = get_reaction_data(tmodel, solution.raw)
# Create the flux dict
flux_dict = this_flux_solution[[i for i in
                                kmodel.reactions.keys()]].to_dict()

# Create a concentration dict with consistent names
variable_names = tmodel.log_concentration.list_attr('name')
metabolite_ids = tmodel.log_concentration.list_attr('id')

temp_concentration_dict = np.exp(solution.raw[variable_names]).to_dict()

# Map concentration names
mapping_dict = {k: sanitize_cobra_vars(v) for k, v in zip(variable_names,
                                                          metabolite_ids)}
concentration_dict = {mapping_dict[k]: v for k, v in
                      temp_concentration_dict.items()}



"""
Sample the kinetic parameters based on linear stability
"""
kmodel.prepare(mca=True)
# Compile MCA functions
kmodel.compile_mca(sim_type=QSSA)

# Initialize parameter sampler
sampling_parameters = SimpleParameterSampler.Parameters(n_samples=10)
sampler = SimpleParameterSampler(sampling_parameters)

# Sample the model
parameter_population = sampler.sample(kmodel, flux_dict,
                                      concentration_dict)
parameter_population = ParameterValuePopulation(parameter_population, kmodel)

# Perform resampling
resampler = SimpleResampler(sampling_parameters)
parameters_to_resample = [kmodel.parameters.km_substrate_ENO, ]
resampled_population = resampler.sample(kmodel, flux_dict,
                                        concentration_dict, parameters_to_resample,
                                        parameter_population._data)
# ToDo user theParameterValuePopulation and ParameterValues objects as input and output of the sampling functions
resampled_population = ParameterValuePopulation(resampled_population, kmodel)

# Calculate the output for the two parameter populations

