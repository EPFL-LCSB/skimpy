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
This tutorial demonstrates modal analysis, MCA, large parameter changes, and basins of attraction
at the example of a metabolic model
"""
from pytfa.io.json import load_json_model
from skimpy.io.yaml import  load_yaml_model
from skimpy.analysis.oracle.load_pytfa_solution import load_concentrations, load_fluxes
from skimpy.analysis.ode.utils import make_flux_fun
from skimpy.analysis.modal import modal_matrix
from skimpy.core.parameters import ParameterValues
from skimpy.utils.namespace import *
from skimpy.utils.tabdict import TabDict

from skimpy.viz.modal import plot_modal_matrix
from skimpy.viz.controll_coefficients import plot_control_coefficients
from skimpy.viz.plotting import timetrace_plot

# Units of the parameters are muM and hr
CONCENTRATION_SCALING = 1e6
TIME_SCALING = 1 # 1hr to 1min
DENSITY = 1200 # g/L
GDW_GWW_RATIO = 0.3 # Assumes 70% Water

kmodel =  load_yaml_model('./../../models/varma_strain_1.yml')
tmodel = load_json_model('./../../models/tfa_varma.json')

# Reference steady-state data
ref_solution = pd.read_csv('./../../data/tfa_reference_strains.csv',
                           index_col=0).loc['strain_1',:]

ref_concentrations = load_concentrations(ref_solution, tmodel, kmodel,
                                         concentration_scaling=CONCENTRATION_SCALING)
ref_fluxes = load_fluxes(ref_solution, tmodel, kmodel,
                               density=DENSITY,
                               ratio_gdw_gww=GDW_GWW_RATIO,
                               concentration_scaling=CONCENTRATION_SCALING,
                               time_scaling=TIME_SCALING)

parameter_values = {p.symbol:p.value for p in kmodel.parameters.values()}
parameter_values = ParameterValues(parameter_values, kmodel)

"""
Modal analysis
"""
kmodel.prepare()
kmodel.compile_jacobian(sim_type=QSSA,ncpu=8)
M = modal_matrix(kmodel,ref_concentrations,parameter_values)

plot_modal_matrix(M,filename='./output/modal_matrix.html',
                  plot_width=800, plot_height=600,
                  clustered=True,
                  backend='svg',
                  )

"""
Metabolic control analysis
"""

# Compile mca with parameter elasticities with respect to Vmaxes
parameter_list = TabDict([(k, p.symbol) for k, p in
                          kmodel.parameters.items() if
                          p.name.startswith('vmax_forward')])

kmodel.compile_mca(sim_type=QSSA,ncpu=8, parameter_list=parameter_list)

flux_control_coeff = kmodel.flux_control_fun(ref_fluxes,
                                             ref_concentrations,
                                             [parameter_values, ])

lac_control_coeff = flux_control_coeff.slice_by('sample',0).loc['LDH_D', :]

lac_control_coeff.index = [v.replace('vmax_forward_','') for v in lac_control_coeff.index ]
plot_control_coefficients(lac_control_coeff,
                          filename='./output/lac_control_coeff.html',
                          backend='svg',
                          )

"""
Large parameter perturbations 
"""

kmodel.compile_ode(sim_type=QSSA,ncpu=8)
# make function to calculate the fluxes
flux_fun = make_flux_fun(kmodel, QSSA)

for k in kmodel.initial_conditions:
    kmodel.initial_conditions[k] = ref_concentrations[k]

desings = {'vmax_forward_LDH_D': 2.0,
           'vmax_forward_GAPD': 2.0}

fluxes = []
for p,v in desings.items():
    kmodel.parameters = parameter_values
    # Implement parameter desing
    kmodel.parameters[p].value = kmodel.parameters[p].value*v

    sol = kmodel.solve_ode(np.logspace(-9,0, 1000),
                             solver_type='cvode')

    # Calculate fluxes:
    this_fluxes = []
    for i, concentrations in sol.concentrations.iterrows():
        t_fluxes = flux_fun(concentrations, parameters=parameter_values)
        this_fluxes.append(t_fluxes)

    # Make it a DataFrame
    this_fluxes = pd.DataFrame(this_fluxes)/ref_fluxes
    this_fluxes.index = sol.time
    fluxes.append(this_fluxes)

ldh_fluxes = np.array([fluxes[0]['LDH_D'].values, fluxes[1]['LDH_D'].values ]).T

timetrace_plot(sol.time,ldh_fluxes ,
               filename='output/ldh_flux.html',
               legend=['2 x [LDH]','2 x [GAPD]'],
               backend='svg'
               )


"""
Basins of attraction
"""

