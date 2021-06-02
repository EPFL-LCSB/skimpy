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
from pytfa.io.json import load_json_model

import numpy as np
import pandas as pd
from skimpy.core.modifiers import *
from skimpy.io.yaml import load_yaml_model
from skimpy.core.reactor import Reactor
from skimpy.analysis.oracle.load_pytfa_solution import load_concentrations, load_fluxes
from skimpy.viz.plotting import timetrace_plot

from pytfa.io.json import load_json_model
from skimpy.io.yaml import load_yaml_model
from skimpy.analysis.oracle.load_pytfa_solution import load_concentrations
from skimpy.core.parameters import load_parameter_population
from skimpy.simulations.reactor import make_batch_reactor
from skimpy.core.solution import ODESolutionPopulation
from skimpy.utils.namespace import *

from skimpy.viz.escher import animate_fluxes, plot_fluxes


import numpy as np
import pandas as pd

WITH_ANIMATION = False # This can take some time

"""
Set up batch reactor
"""

reactor = make_batch_reactor('single_species.yaml')
reactor.compile_ode(add_dilution=False)

"""
"""
path_to_kmodel = './../../models/kin_varma.yml'
path_to_tmodel = './../../models/tfa_varma.json'

# load models
tmodel = load_json_model(path_to_tmodel)
kmodel = load_yaml_model(path_to_kmodel)


reference_solutions = pd.read_csv('./../../data/tfa_reference_strains.csv', index_col=0)
ref_concentrations = load_concentrations(reference_solutions.loc['strain_1'], tmodel, kmodel,
                                                  concentration_scaling=reactor.concentration_scaling)

reactor.initial_conditions['biomass_strain_1'] = 0.1e12 # Number of cells

reactor.initialize(ref_concentrations, 'strain_1')

"""
Solve 
"""

sol = reactor.solve_ode(np.linspace(0, 10.0, 1000),
                        solver_type='cvode',
                        rtol=1e-9,
                        atol=1e-9,
                        max_steps=1e9,
                        )

"""
Plot results 
"""

species = [s for s in sol.concentrations.columns if not s in ['biomass_strain_1',]]
timetrace_plot(sol.time, sol.concentrations[species].values/reactor.concentration_scaling,
               filename='output_single/time_response.html',
               legend=species,
               x_label='time [h]',
               y_label='concentrations [M]',
               backend='svg',)

MASS_PER_CELL = 1e-12 #[g/cell]
species = ['biomass_strain_1', ]
timetrace_plot(sol.time, sol.concentrations[species].values*MASS_PER_CELL,
               filename='output_single/time_response_biomass.html',
               legend=species,
               x_label='time [h]',
               y_label='biomass [g]',
               legend_location='top_left',
               backend='svg',)

# Calc growth rate
species = ['biomass_strain_1', ]
delta_x =  (sol.concentrations[species].values[1:,:] -  sol.concentrations[species].values[:-1,:])
delta_t = np.array(  [(sol.time[1:] -  sol.time[:-1])]).T
x = sol.concentrations[species].values[:-1,:]
mu = delta_x/delta_t/x
timetrace_plot(sol.time[:-1], mu,
               filename='output_single/time_response_growth.html',
               legend=species,
               x_label='time [h]',
               y_label='growth rate [1/h]',
               backend='svg',
               )

# Medium
species = list(reactor.medium.keys())
timetrace_plot(sol.time, sol.concentrations[species].values/reactor.concentration_scaling,
               filename='output_single/time_response_medium.html',
               legend=species,
               x_label='time [h]',
               y_label='concentrations [M]',
               backend='svg',)

"""
Compute fluxes Strain 1
"""

from skimpy.analysis.ode.utils import make_flux_fun

fluxes_strain_1 = make_flux_fun(reactor.models.strain_1, QSSA)
fluxes = []
model_params = reactor.models.strain_1.parameters
parameters = {str(p.symbol):p.value for p in model_params.values() if not p.value is None}
for i, concentrations in sol.concentrations.iterrows():
    this_fluxes = fluxes_strain_1(concentrations, parameters=parameters)
    fluxes.append(this_fluxes)

# Make it a DataFrame
fluxes = pd.DataFrame(fluxes)/reactor.flux_scaling
fluxes.index = sol.time

"""
Animate the fluxes solutions in a escher map
"""
if WITH_ANIMATION:
    animate_fluxes(fluxes[::10],
                   './../../data/varma_map.json',
                   outputfile='output_single/animation_fluxes_strain_1.mp4',
                   height=600,
                   width=800,
                   time_interval_ms=100,
                   min_flux = -10,
                   max_flux = 10,)

plot_fluxes( fluxes.iloc[0],
             './../../data/varma_map.json',
             output_file='output_single/map_1.html',
             height=600,
             width=800,
             min_flux = -10,
             max_flux = 10,)

plot_fluxes( fluxes.iloc[100],
             './../../data/varma_map.json',
             output_file='output_single/map_2.html',
             height=600,
             width=800,
             min_flux = -10,
             max_flux = 10,)
