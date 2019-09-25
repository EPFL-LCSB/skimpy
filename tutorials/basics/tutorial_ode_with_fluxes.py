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

# Test models
import numpy as np
from skimpy.core import *
from skimpy.mechanisms import *
from skimpy.analysis.ode.utils import make_flux_fun
from skimpy.utils.namespace import *

name = 'pfk'
metabolites = ReversibleMichaelisMenten.Reactants(substrate = 'A',
                                                   product = 'B')

## QSSA Method
parameters = ReversibleMichaelisMenten.Parameters(
    vmax_forward = 1.0,
    k_equilibrium = 1.5,
    km_substrate = 10.0,
    km_product = 10.0,
    total_enzyme_concentration = 1.0,
)

pfk = Reaction(name=name,
               mechanism = ReversibleMichaelisMenten,
               reactants=metabolites,
               )

this_model = KineticModel()
this_model.add_reaction(pfk)
this_model.parametrize_by_reaction({pfk.name:parameters})


## Elementary rate method
this_model.compile_ode(sim_type = ELEMENTARY)

this_model.initial_conditions['A'] = 10.0
this_model.initial_conditions['B']= 1.0
this_model.initial_conditions['pfk'] = 1.0

this_sol_full = this_model.solve_ode(np.linspace(0.0, 100.0, 1000), solver_type='cvode')

calc_fluxes = make_flux_fun(this_model, ELEMENTARY)

steady_state_fluxes = calc_fluxes(this_sol_full.concentrations.iloc[-1])

