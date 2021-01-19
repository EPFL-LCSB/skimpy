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



"""
This Tutorial walks you through the ORACLE workflow and show case several analysis methods
implemented in SkiMpy
Each section contains a short introduction on the purpose and pitfalls 
Further the key parameters of the step are shown in the beginning of each step
"""


"""
pyTFA sampling
"""
NUM_TFA_SAMPLES = 100

# Note: This a preprocessed model curated to exhibit a single
# flux directionality profile (FDP) and is stripped of all integer variables
# For details on enumeration of FDPs and sampling pre-processing please refer to
# pytfa/tutorials/ XX and XX

path_to_tmodel = './../models/tfa_varma.json'
tmodel = load_json_model(path_to_tmodel)

GLPK= 'optlang-glpk'
CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
tmodel.solver = GLPK

tmodel.solver.configuration.tolerances.feasibility = 1e-9
tmodel.solver.configuration.tolerances.optimality  = 1e-9
tmodel.solver.configuration.tolerances.integrality = 1e-9

# Test feasiblity
print(tmodel.optimize())

samples = sample(continuous_model, NUM_TFA_SAMPLES, method='achr')
samples.to_csv('./output/samples.csv'.format())

"""
Parameter sampling and pruning
"""




"""
Modal analysis
"""



"""
Basins of attraction
"""



"""
Metabolic control analysis 
"""







