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


import numpy as np
import pandas as pd

from skimpy.core import *
from skimpy.mechanisms import *
from skimpy.viz.plotting import timetrace_plot

Degradation = make_irrev_massaction([-1])

"""
Here we model a GFP stripe forming genetic gene circuit based scheme on the SI in Schaerli et al, 2014
Schaerli, Y., Munteanu, A., Gili, M. et al. A unified design space of synthetic stripe-forming networks. 
Nat Commun 5, 4905 (2014). https://doi.org/10.1038/ncomms5905
"""

k_deg = 4 # This param is guessed until a stripe was observed not found in the SI of above

#  Expression of SP6
reactants = SimpleRegulatedGeneExpression.Reactants(product='SP6',)
regulators = SimpleRegulatedGeneExpression.SingleRegulator(factor_1='arabinose',)
r1 = Reaction(name='SP6_expression',
              mechanism=SimpleRegulatedGeneExpression,
              inhibitors=regulators,
              reactants=reactants,)
p_r1 = SimpleRegulatedGeneExpression.Parameters(basal_transcription_rate=0.862,
                                                transcription_rate_factor_1=29.7,
                                                hill_coefficient_factor_1=1,
                                                kd_factor_1=1/111.0,)
# Degradation of TetR
reactants = Degradation.Reactants(substrate1='SP6',)
r2 = Reaction(name='SP6_degradation',
              mechanism=Degradation,
              reactants=reactants,)
p_r2 = Degradation.Parameters(vmax_forward=k_deg,)

# Expression of Lacl
reactants = SimpleRegulatedGeneExpression.Reactants(product='lacl',)
regulators = SimpleRegulatedGeneExpression.SingleRegulator(factor_1='SP6',)
r3 = Reaction(name='lacl_expression',
              mechanism=SimpleRegulatedGeneExpression,
              inhibitors=regulators,
              reactants=reactants,)
p_r3 = SimpleRegulatedGeneExpression.Parameters(basal_transcription_rate=4.76,
                                                transcription_rate_factor_1=48.5,
                                                hill_coefficient_factor_1=1,
                                                kd_factor_1=1/0.0645,)
# Degradation of TetR
reactants = Degradation.Reactants(substrate1='lacl',)
r4 = Reaction(name='lacl_degradation',
              mechanism=Degradation,
              reactants=reactants,)
p_r4 = Degradation.Parameters(vmax_forward=k_deg,)


# Expression of GFP
reactants = SimpleRegulatedGeneExpression.Reactants(product='GFP',)
regulators = SimpleRegulatedGeneExpression.Regulators(factor_1='SP6', factor_2='lacl')
r5 = Reaction(name='GFP_expression',
              mechanism=SimpleRegulatedGeneExpression,
              inhibitors=regulators,
              reactants=reactants,)
p_r5 = SimpleRegulatedGeneExpression.Parameters(basal_transcription_rate=0.0,
                                                transcription_rate_factor_1=1.78e3,
                                                hill_coefficient_factor_1=1.4, # This is supposed to 14 in Sup mat from above
                                                kd_factor_1=1.0,
                                                transcription_rate_factor_2=0,
                                                hill_coefficient_factor_2=2.0,
                                                kd_factor_2=1/0.33,
                                                transcription_rate_factor_1_factor_2=72.4,
                                                cooperativity=46.5,
                                                )
# Degradation of GFP
reactants = Degradation.Reactants(substrate1='GFP',)
r6 = Reaction(name='GFP_degradation',
              mechanism=Degradation,
              reactants=reactants,)
p_r6 = Degradation.Parameters(vmax_forward=k_deg,)


kmodel = KineticModel()
for r in [r1,r2,r3,r4,r5,r6]:
    kmodel.add_reaction(r)
kmodel.parametrize_by_reaction({r1.name:p_r1, r2.name:p_r2, r3.name:p_r3,
                                    r4.name:p_r4, r5.name:p_r5, r6.name:p_r6, })

kmodel.repair()
kmodel.compile_ode(sim_type = QSSA)

kmodel.initial_conditions['arabinose'] = 10.

sol = kmodel.solve_ode(np.linspace(0.0, 100.0, 1000), solver_type='cvode')
sol.plot('output/time_response.html', backend='svg')


input = np.logspace(-5,1,100)
output = []
for ara in input:
    kmodel.initial_conditions['arabinose'] = ara
    sol = kmodel.solve_ode(np.linspace(0.0, 1000.0, 1000), solver_type='cvode')
    # Get last point
    output.append(sol.concentrations.loc[999,['SP6','lacl','GFP',]])

output = np.array(output)

timetrace_plot(input, output,
               filename='output/input_output.html',
               legend=['SP6','lacl','GFP',],
               x_axis_type="log",
               backend='svg',
               )

