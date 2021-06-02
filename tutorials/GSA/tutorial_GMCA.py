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
from pytfa.io.viz import get_reaction_data, export_reactions_for_escher

from skimpy.core import *
from skimpy.mechanisms import *
from skimpy.utils.namespace import *
from skimpy.sampling.simple_resampler import SimpleResampler
from skimpy.sampling.simple_parameter_sampler import SimpleParameterSampler
from skimpy.core.parameters import ParameterValuePopulation
from skimpy.io.generate_from_pytfa import FromPyTFA
from skimpy.utils.general import sanitize_cobra_vars
from skimpy.utils.tabdict import TabDict

from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, \
    load_concentrations, load_equilibrium_constants

from skimpy.analysis.oracle import *
from skimpy.viz.sensitivty import plot_sobol_coefficients

"""
Parameters
"""
CONCENTRATION_SCALING = 1e6 # 1 mol to 1 mmol
TIME_SCALING = 1 # 1hr
# Parameters of the E. Coli cell
DENSITY = 1200 # g/L
GDW_GWW_RATIO = 0.3 # Assumes 70% Water

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
# with minium fluxes of 1.0
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
               'GLYCK':       (1.0, 100.0),
               'GLYCK2':      (1.0, 100.0),
               'PGK':         (1.0, 100.0),
               'PGM':         (1.0, 100.0),
               'TRSARr':      (2.0, 10.0),
               'Trp_adp':     (-100.0, 100.0),
               'Trp_atp':     (-100.0, 100.0),
               'Trp_h':       (-100.0, 100.0),
               'Trp_h2o':     (-100.0, 100.0),
               'Trp_nad':     (-100.0, 100.0),
               'Trp_nadh':    (-100.0, 100.0)}

for b, (lb, ub) in this_bounds.items():
    tmodel.reactions.get_by_id(b).bounds = (lb, ub)


# Find a solution for this FDP
solution = tmodel.optimize()

# Force a minimal thermodynamic displacement
min_log_displacement = 1e-1
add_min_log_displacement(tmodel, min_log_displacement)

# Find a solution for the model
solution = tmodel.optimize()

# Export solution for Escher
export_reactions_for_escher(tmodel,solution.raw, 'output/reference_fluxes.csv')

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

fluxes = load_fluxes(solution.raw, tmodel, kmodel,
                     density=DENSITY,
                     ratio_gdw_gww=GDW_GWW_RATIO,
                     concentration_scaling=CONCENTRATION_SCALING,
                     time_scaling=TIME_SCALING)

concentrations = load_concentrations(solution.raw, tmodel, kmodel,
                                     concentration_scaling=CONCENTRATION_SCALING)

# Fetch equilibrium constants
load_equilibrium_constants(solution.raw, tmodel, kmodel,
                           concentration_scaling=CONCENTRATION_SCALING,
                           in_place=True)

"""
Sample the kinetic parameters based on linear stability
"""
kmodel.prepare(mca=True)

# Compile MCA functions
parameter_list = TabDict([(k, p.symbol) for k, p in
                          kmodel.parameters.items() if
                          p.name.startswith('vmax_forward')])

kmodel.compile_mca(sim_type=QSSA, parameter_list=parameter_list)

# Define a function that calculates model output (flux control coefficients)
# from sampled model parameters

# Choose a subset of Vmax's to analyse
outputs_to_analyse = ['vmax_forward_ENO','vmax_forward_PGK']

def calculate_model_output(model_parameters):
    # calculate all flux control coefficients
    flux_control_coeff_0 = kmodel.flux_control_fun(fluxes,
                                                concentrations,
                                                model_parameters._data)

    # Choose a particular reaction we want to analyse
    df = flux_control_coeff_0.slice_by('flux','ENO')
    # Select only some of the Vmax's
    df = df[df.index.isin(outputs_to_analyse)]
    return df

# Initialize parameter sampler
sampling_parameters = SimpleParameterSampler.Parameters(n_samples=200)
sampler = SimpleParameterSampler(sampling_parameters)

# Sample the model (A matrix)
parameter_population = sampler.sample(kmodel, fluxes, concentrations)
parameter_population = ParameterValuePopulation(parameter_population, kmodel)

# Construct B,C matrices given parameter_population and parameters_to_resample
parameters_to_resample = [[kmodel.parameters.km_substrate_ENO, ],
                          [kmodel.parameters.km_product_ENO, ],
                          [kmodel.parameters.km_substrate_PGM, ],
                          [kmodel.parameters.km_product_PGM, ],
                          [kmodel.parameters.km_substrate1_PGK, ],
                          [kmodel.parameters.km_substrate2_PGK, ],
                          [kmodel.parameters.km_product1_PGK, ],
                          [kmodel.parameters.km_product2_PGK, ]]

df_si = pd.DataFrame(index=outputs_to_analyse)
df_st = pd.DataFrame(index=outputs_to_analyse)

for this_parameter in parameters_to_resample:
    outer_set = [i for i in kmodel.parameters.values() if i not in this_parameter]

    # Resample the initial parameter_population according to this_parameter
    resampler = SimpleResampler(sampling_parameters)

    resampled_population_B = resampler.sample(kmodel, fluxes, concentrations,
                                              outer_set,
                                              parameter_population._data)
    resampled_population_B = ParameterValuePopulation(resampled_population_B, kmodel)

    resampled_population_C = resampler.sample(kmodel, fluxes, concentrations,
                                              this_parameter,
                                              parameter_population._data)
    resampled_population_C = ParameterValuePopulation(resampled_population_C, kmodel)


    # Calculate output
    df_A = calculate_model_output(parameter_population)
    df_B = calculate_model_output(resampled_population_B)
    df_C = calculate_model_output(resampled_population_C)

    # Calculate sensitivity indices for each flux control coefficient/model output
    df_sensitivity_indices = pd.DataFrame(columns=['si', 'st'])

    for index in df_A.index:
        n = len(df_A.loc[index])
        f0 = np.sum(df_A.loc[index])/n

        si = (np.dot(df_A.loc[index], df_B.loc[index])/n - f0**2) / (np.dot(df_A.loc[index], df_A.loc[index])/n - f0**2)
        st = (np.dot(df_A.loc[index], df_C.loc[index])/n - f0**2) / (np.dot(df_A.loc[index], df_A.loc[index])/n - f0**2)

        df_sensitivity_indices.loc[index] = pd.Series({'si':si, 'st':st})

    df_si[this_parameter[0].symbol] = df_sensitivity_indices['si']
    df_st[this_parameter[0].symbol] = df_sensitivity_indices['st']

# TODO: concatenate A,B,C in to one parameter pop to have only one function call
# to calculate_model_output
# TODO: specify parameters_to_resample as a dict, so that we can resample groups
# of parameters and let the user name these groups

# Plot
plot_sobol_coefficients(df_si.loc['vmax_forward_ENO'],
                        df_st.loc['vmax_forward_ENO'],
                        filename='output/eno_si_st.html',
                        colors=['#3399FF','#205f99'],
                        backend='svg')
