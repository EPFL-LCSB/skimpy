# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2021 Laboratory of Computational Systems Biotechnology (LCSB),
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

"""
This Tutorial illustrates how signaling models can be build using skimpy
"""

import numpy as np
import pandas as pd

from skimpy.core import *
from skimpy.mechanisms import *
from skimpy.core.modifiers import InhibitionModifier
from skimpy.viz.plotting import timetrace_plot

# Get Class for the MichaelisMenten mechanism
MichaelisMenten= make_irrev_m_n_michaelis_menten([-1, 1])

# Receptor
metabolites = MichaelisMenten.Reactants(substrate1='MAPKKK',  product1='MAPKKKP',)
r1 = Reaction(name='Receptor', mechanism=MichaelisMenten, reactants=metabolites,)
# With
r1.modifiers['MAPKP'] = InhibitionModifier('MAPKP', k_inhibition=0.1)

p_r1 = MichaelisMenten.Parameters(vmax_forward=50, km_substrate1=10.0,)

# Phosphatases 1
metabolites = MichaelisMenten.Reactants(substrate1='MAPKKKP',  product1='MAPKKK',)
r2 = Reaction(name='Phosphatases_1', mechanism=MichaelisMenten, reactants=metabolites,)
p_r2 = MichaelisMenten.Parameters(vmax_forward=0.25, km_substrate1=8.0,)

# MAPKKK-P
metabolites = MichaelisMenten.Reactants(substrate1='MAPKK',  product1='MAPKKP',)
r3 = Reaction(name='MAPKKKP', mechanism=MichaelisMenten, reactants=metabolites, enzyme='MAPKKKP')
p_r3 = MichaelisMenten.Parameters(kcat_forward=0.025, km_substrate1=15.0,)

# Phosphatases 2
metabolites = MichaelisMenten.Reactants(substrate1='MAPKKP',  product1='MAPKK',)
r4 = Reaction(name='Phosphatases_2', mechanism=MichaelisMenten, reactants=metabolites,)
p_r4 = MichaelisMenten.Parameters(vmax_forward=0.75, km_substrate1=15.0,)

# MAPKK-P
metabolites = MichaelisMenten.Reactants(substrate1='MAPK',  product1='MAPKP',)
r5 = Reaction(name='MAPKKP', mechanism=MichaelisMenten, reactants=metabolites, enzyme='MAPKKP')
p_r5 = MichaelisMenten.Parameters(kcat_forward=0.025, km_substrate1=15.0,)

# Phosphatases 3
metabolites = MichaelisMenten.Reactants(substrate1='MAPKP',  product1='MAPK',)
r6 = Reaction(name='Phosphatases_3', mechanism=MichaelisMenten, reactants=metabolites,)
p_r6 = MichaelisMenten.Parameters(vmax_forward=0.5, km_substrate1=15.0,)

this_model = KineticModel()
for r in [r1,r2,r3,r4,r5,r6]:
    this_model.add_reaction(r)
this_model.parametrize_by_reaction({r1.name:p_r1, r2.name:p_r2, r3.name:p_r3,
                                    r4.name:p_r4, r5.name:p_r5, r6.name:p_r6, })
this_model.compile_ode(sim_type = QSSA)

this_model.initial_conditions['MAPKKK'] = 100.0
this_model.initial_conditions['MAPKK'] = 300.0
this_model.initial_conditions['MAPK'] = 300.0

this_sol_full = this_model.solve_ode(np.linspace(0.0, 10000.0, 1000),
                                     solver_type='cvode')

this_sol_full.plot('output/time_signaling_response_oscillations.html',
                   backend='svg',
                   )







