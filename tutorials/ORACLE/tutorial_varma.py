# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2020 Laboratory of Computational Systems Biotechnology (LCSB),
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

import pandas as pd
import numpy as np

"""
This Tutorial walks you through the ORACLE workflow and show case several analysis methods
implemented in SkiMpy
Each section contains a short introduction on the purpose and pitfalls 
Further the key parameters of the step are shown in the beginning of each step
"""


"""
pyTFA sampling
"""

from skimpy.analysis.oracle.minimum_fluxes import MinFLux, \
    MinFLuxVariable

from pytfa.io.json import load_json_model, save_json_model
from pytfa.analysis.sampling import sample

NUM_TFA_SAMPLES = 10

# Note: This a preprocessed model curated to exhibit a single
# flux directionality profile (FDP) and is stripped of all integer variables
# For details on enumeration of FDPs and sampling pre-processing please refer to
# pytfa/tutorials/ XX and XX

path_to_tmodel = './../../models/tfa_varma.json'
tmodel = load_json_model(path_to_tmodel)

GLPK= 'optlang-glpk'
CPLEX = 'optlang-cplex'
tmodel.solver = GLPK

tmodel.solver.configuration.tolerances.feasibility = 1e-9
#tmodel.solver.configuration.tolerances.optimality  = 1e-9

# Define the reactor medium
#Glucose  to 10 g/L = 0.056 mol/l 1/2 Reactor concentration
GLUCOSE_CONCENTRATION = 0.056
PHOSPHATE = 10e-3
CARBON_DIOXIDE = 1e-7
OXYGEN = 8e-3*0.062 # 8 mg/L 1g = 0.062 mol

tmodel.log_concentration.get_by_id('glc-D_e').variable.ub = np.log(GLUCOSE_CONCENTRATION*1.2)
tmodel.log_concentration.get_by_id('glc-D_e').variable.lb = np.log(GLUCOSE_CONCENTRATION*0.8)

tmodel.log_concentration.get_by_id('pi_e').variable.lb = np.log(PHOSPHATE*0.8)
tmodel.log_concentration.get_by_id('pi_e').variable.ub = np.log(PHOSPHATE*1.2)

tmodel.log_concentration.get_by_id('co2_e').variable.ub = np.log(CARBON_DIOXIDE*1.2)
tmodel.log_concentration.get_by_id('co2_e').variable.lb = np.log(CARBON_DIOXIDE*0.8)

tmodel.log_concentration.get_by_id('o2_e').variable.ub = np.log(OXYGEN*1.2)
tmodel.log_concentration.get_by_id('o2_e').variable.lb = np.log(OXYGEN*0.8)

# Enforce glucose transporter displacements
tmodel.thermo_displacement.GLCptspp.variable.ub = -2.0

# Reduce the other
# Constraint non-medium concentrations to be lower then muM
LC_SECRETION = np.log(1e-6)
secretions = [r for r in tmodel.boundary if r.upper_bound <= 0]
for sec in secretions:
    for met in sec.metabolites:
        if met.id in ['h2o_2', 'h_e']:
            continue
        try:
            tmodel.log_concentration.get_by_id(met.id).variable.upper_bound = LC_SECRETION
        except KeyError:
            pass


# Test feasiblity
print(tmodel.optimize())

samples = sample(tmodel, NUM_TFA_SAMPLES, method='achr')
samples.to_csv('./output/samples_fdp1_1000.csv'.format())

"""
Parameter sampling 
"""
# Note this kinetic model was created from a draft model using the pytfa exporter
# See tutorials/io/import_pytfa.py or tutorials/ORACLE/tutorial_toy for reference

from skimpy.io.yaml import load_yaml_model
from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, \
    load_concentrations, load_equilibrium_constants
from skimpy.core.parameters import ParameterValuePopulation, load_parameter_population
from skimpy.sampling.simple_parameter_sampler import SimpleParameterSampler
from skimpy.utils.general import get_stoichiometry


NCPU = 8
N_PARAMETER_SAMPLES = 10
CONCENTRATION_SCALING = 1e6 # 1 mol to 1 mmol
TIME_SCALING = 1 # 1hr
# Parameters of the E. Coli cell
DENSITY = 1200 # g/L
GDW_GWW_RATIO = 0.3 # Assumes 70% Water

flux_scaling_factor = 1e-3 / (GDW_GWW_RATIO / DENSITY) \
                      * CONCENTRATION_SCALING \
                      / TIME_SCALING


path_to_kmodel = './../../models/kin_varma.yml'
path_for_output = './output/paramter_pop_{}.h5'
kmodel = load_yaml_model(path_to_kmodel)

# Perp and compile to sample parameters
kmodel.prepare()
kmodel.compile_jacobian(ncpu=NCPU)

# Constraint glucose transporter KMs
kmodel.reactions.GLCptspp.parameters.km_substrate2.bounds = (14.99, 15.01)




sampler_params = SimpleParameterSampler.Parameters(n_samples=N_PARAMETER_SAMPLES)
sampler = SimpleParameterSampler(sampler_params)

lambda_max_all = []
lambda_min_all = []

S = get_stoichiometry(kmodel, kmodel.reactants).todense()

from skimpy.utils.namespace import *
from skimpy.analysis.ode.utils import make_flux_fun
fluxfun= make_flux_fun(kmodel, QSSA)
fluxes = []

for i, sample in samples.iterrows():
    # Load fluxes and concentrations
    fluxes = load_fluxes(sample, tmodel, kmodel,
                         density=DENSITY,
                         ratio_gdw_gww=GDW_GWW_RATIO,
                         concentration_scaling=CONCENTRATION_SCALING,
                         time_scaling=TIME_SCALING)
    concentrations = load_concentrations(sample, tmodel, kmodel,
                                         concentration_scaling=CONCENTRATION_SCALING)

    # Calibrate the glucose transporter to have a Vmax of 10mmol/gDW/h
    pep_c = concentrations['pep_c']
    pyr_c = concentrations['pyr_c']
    g6p_c = concentrations['g6p_c']

    kmodel.reactions.GLCptspp.parameters.km_substrate1.bounds = (pep_c*1e-4, pep_c*1e-3)
    kmodel.reactions.GLCptspp.parameters.km_product1.bounds = (pyr_c*1.0, pyr_c*1e2)
    kmodel.reactions.GLCptspp.parameters.km_product2.bounds = (g6p_c*1.0, g6p_c*1e2)

    # Unsaturated Accetate transporter
    ac_c = concentrations['ac_c']
    kmodel.reactions.ACt2rpp.parameters.km_product.bounds = (ac_c*1e2, ac_c*1e3)

    # ATP maintainance to be saturated
    atp_c = concentrations['atp_c']
    kmodel.reactions.ATPM.parameters.km_substrate1.bounds = (atp_c*1e-4, atp_c*1e-3)

    # ATPS to be stuarated in exp_phase
    pi_c = concentrations['pi_c']
    adp_c = concentrations['adp_c']
    kmodel.reactions.ATPS4rpp.parameters.km_substrate1.bounds = (pi_c*1e-4, pi_c*1e-3)
    kmodel.reactions.ATPS4rpp.parameters.km_substrate2.bounds = (adp_c*1e-4, adp_c*1e-3)

    # Fetch equilibrium constants
    load_equilibrium_constants(sample, tmodel, kmodel,
                               concentration_scaling=CONCENTRATION_SCALING,
                               in_place=True)


    # Test weather the model stoichiometry is consistent with the steady-state
    dxdt = S.dot(fluxes)
    if np.any(abs(dxdt) > 1e-9*flux_scaling_factor):
        raise RuntimeError('dxdt for idx {} not equal to 0'
                           .format(np.where(abs(dxdt) > 1e-9*flux_scaling_factor)))

    indep_reactants = [kmodel.reactants.iloc(i)[0] for i in kmodel.independent_variables_ix]
    dxdt_red = kmodel.reduced_stoichiometry.todense().dot(fluxes)
    if np.any(abs(dxdt_red) > 1e-9*flux_scaling_factor):
        raise RuntimeError('dxdt reduced for idx {} not equal to 0'
                           .format(np.where(abs(dxdt_red) > 1e-9*flux_scaling_factor)))

    # Generate sampels and fetch slowest and fastest eigenvalues
    params, lamda_max, lamda_min = sampler.sample(kmodel, fluxes, concentrations,
                                                  only_stable=True,
                                                  min_max_eigenvalues=True)
    lambda_max_all.append(pd.DataFrame(lamda_max))
    lambda_min_all.append(pd.DataFrame(lamda_min))

    params_population = ParameterValuePopulation(params, kmodel)
    params_population.save(path_for_output.format(i))

    # Test if the resulting sets are NV=0
    fluxes_1 = fluxfun(concentrations, parameters=params_population['0'])
    fluxes_1 = pd.Series(fluxes_1)
    dxdt = S.dot(fluxes_1[kmodel.reactions])
    if np.any(abs(dxdt) > 1e-9*flux_scaling_factor):
        raise RuntimeError('dxdt for idx {} not equal to 0'
                           .format(np.where(abs(dxdt) > 1e-9*flux_scaling_factor)))



# Process df and save dataframe
lambda_max_all = pd.concat(lambda_max_all, axis=1)
lambda_min_all = pd.concat(lambda_min_all, axis=1)


"""
Prune parameters based on the time scales
"""

MAX_EIGENVALUES = -5/(20/60)    # 1/hr
MIN_EIGENVALUES = -3e11

# Prune parameter based on eigenvalues
is_selected = (lambda_max_all < MAX_EIGENVALUES ) & (lambda_min_all > MIN_EIGENVALUES )
is_selected.columns = range(lambda_max_all.shape[1])

fast_parameters = []
fast_index = []

for i, row in is_selected.T.iterrows():
    if any(row):
        fast_models = np.where(np.array(row))[0]
        # Load the respective solutions
        parameter_population = load_parameter_population(path_for_output.format(i))
        fast_parameters.extend([parameter_population._data[k] for k in fast_models])
        fast_index.extend(["{},{}".format(i,k) for k in fast_models])

# Generate a parameter population file
parameter_population = ParameterValuePopulation(fast_parameters,
                                           kmodel=kmodel,
                                           index=fast_index)
parameter_population.save('pruned_parameters.hdf5')



    