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
from skimpy.sampling.simple_parameter_sampler import SimpleParameterSampler

from skimpy.core.parameters import ParameterValuePopulation

from skimpy.core.solution import ODESolutionPopulation
from skimpy.io.generate_from_pytfa import FromPyTFA
from skimpy.utils.general import sanitize_cobra_vars
from skimpy.utils.tabdict import TabDict
from skimpy.analysis.ode.utils import make_flux_fun

from skimpy.analysis.oracle import *

from optlang.exceptions import SolverError

CONCENTRATION_SCALING = 1e6 # 1 mol to 1 mmol
TIME_SCALING = 1 # 1hr
# Parameters of the E. Coli cell
DENSITY = 1200 # g/L
GDW_GWW_RATIO = 0.3 # Assumes 70% Water

flux_scaling_factor = 1e-3 * (GDW_GWW_RATIO * DENSITY) \
                      * CONCENTRATION_SCALING \
                      / TIME_SCALING

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
this_pytfa_model = pytfa.ThermoModel(thermo_data, this_cobra_model)

CPLEX = 'optlang-cplex'
this_pytfa_model.solver = CPLEX

# TFA conversion
this_pytfa_model.prepare()
this_pytfa_model.convert(add_displacement=True)


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
               'GLYCK2':      (1e-3, 100.0),
               'PGK':         (1e-3, 100.0),
               'PGM':         (2.0, 2.0),
               'TRSARr':      (2.0, 10.0),
               'Trp_adp':     (-100.0, 100.0),
               'Trp_atp':     (-100.0, 100.0),
               'Trp_h':       (-100.0, 100.0),
               'Trp_h2o':     (-100.0, 100.0),
               'Trp_nad':     (-100.0, 100.0),
               'Trp_nadh':    (-100.0, 100.0)}

for k,v in this_bounds.items():
    this_pytfa_model.reactions.get_by_id(k).bounds = v

# Find a solution for this FDP
solution = this_pytfa_model.optimize()

# Force a minimal thermodynamic displacement
min_log_displacement = 1e-1
add_min_log_displacement(this_pytfa_model, min_log_displacement)

# Find a solution for the model
solution = this_pytfa_model.optimize()

this_pytfa_model.thermo_displacement.PGM.variable.lb = np.log(1e-1)
this_pytfa_model.thermo_displacement.PGM.variable.ub = np.log(1e-1)

solution = this_pytfa_model.optimize()


"""
Get a Kinetic Model
"""
# Generate the KineticModel

# Define the molecules that should be considered small-molecules
# These molecules will not be accounted explicitly in the kinetic mechanism as
# substrates and products
small_molecules = ['h_c', 'h_e']

model_gen = FromPyTFA(small_molecules=small_molecules)
this_skimpy_model = model_gen.import_model(this_pytfa_model, solution.raw)

"""
Load the reference solutions
"""

fluxes = load_fluxes(solution.raw, this_pytfa_model, this_skimpy_model,
                     density=DENSITY,
                     ratio_gdw_gww=GDW_GWW_RATIO,
                     concentration_scaling=CONCENTRATION_SCALING,
                     time_scaling=TIME_SCALING)

concentrations = load_concentrations(solution.raw, this_pytfa_model, this_skimpy_model,
                                     concentration_scaling=CONCENTRATION_SCALING)

load_equilibrium_constants(solution.raw, this_pytfa_model, this_skimpy_model,
                           concentration_scaling=CONCENTRATION_SCALING,
                           in_place=True)

"""
Sample the kinetic parameters based on linear stablity 
"""
this_skimpy_model.parameters.km_substrate_ENO.bounds = (1e-4, 1e-3)
this_skimpy_model.parameters.km_product_ENO.bounds = (1e-4, 1e-3)

this_skimpy_model.prepare(mca=True)
# Compile MCA functions
this_skimpy_model.compile_mca(sim_type=QSSA)

# Initialize parameter sampler
sampling_parameters = SimpleParameterSampler.Parameters(n_samples=1)
sampler = SimpleParameterSampler(sampling_parameters)

# Sample the model
parameter_population = sampler.sample(this_skimpy_model,
                                      fluxes,
                                      concentrations)
parameter_population = ParameterValuePopulation(parameter_population, this_skimpy_model, index=range(1))

"""
Calculate control coefficients
"""
parameter_list = TabDict([(k, p.symbol) for k, p in
                          this_skimpy_model.parameters.items() if
                          p.name.startswith('vmax_forward')])

this_skimpy_model.compile_mca(mca_type=SPLIT,sim_type=QSSA, parameter_list=parameter_list)

flux_control_coeff = this_skimpy_model.flux_control_fun(fluxes,
                                                        concentrations,
                                                        parameter_population)


"""
Integrate the ODEs
"""

this_skimpy_model.compile_ode(sim_type=QSSA)

this_skimpy_model.initial_conditions = TabDict([(k,v) for k,v in concentrations.iteritems()])


solutions = []

this_parameters = parameter_population[0]
vmax_glyck2 = this_parameters['vmax_forward_GLYCK2']

# For each solution calulate the fluxes
calc_fluxes = make_flux_fun(this_skimpy_model, QSSA)

fluxes_pgm_expression = []

this_skimpy_model.parameters = this_parameters

for rel_e in np.logspace(-3, 3, 100):
    this_parameters['vmax_forward_GLYCK2'] = vmax_glyck2*rel_e

    this_skimpy_model.parameters = this_parameters

    sol = this_skimpy_model.solve_ode(np.linspace(0.0, 1.0, 1000),
                                      solver_type='cvode',
                                      rtol=1e-9,
                                      atol=1e-9,
                                      max_steps=1e9,)

    solutions.append(sol)

    steady_state_fluxes = calc_fluxes(sol.concentrations.iloc[-1], parameters=this_parameters)
    fluxes_pgm_expression.append(steady_state_fluxes)

fluxes_pgm_expression = pd.DataFrame(fluxes_pgm_expression)/flux_scaling_factor
