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

import pytfa
from pytfa.io import import_matlab_model, load_thermoDB
from pytfa.io.viz import get_reaction_data

from skimpy.core import *
from skimpy.mechanisms import *
from skimpy.utils.namespace import *
from skimpy.sampling.simple_parameter_sampler import SimpleParameterSampler
from skimpy.core.solution import ODESolutionPopulation
from skimpy.io.generate_from_pytfa import FromPyTFA
from skimpy.utils.general import sanitize_cobra_vars
from skimpy.utils.tabdict import TabDict

from skimpy.analysis.oracle import *

""" 
Import and curate a model
"""

this_cobra_model = import_matlab_model('../../models/toy_model.mat', 'model')

""" 
Make tfa model
"""

# Convert to a thermodynamics model
thermo_data = load_thermoDB('../../data/thermo_data.thermodb')
this_pytfa_model = pytfa.ThermoModel(thermo_data, this_cobra_model)

GLPK = 'optlang-glpk'
this_pytfa_model.solver = GLPK

## TFA conversion
this_pytfa_model.prepare()
this_pytfa_model.convert(add_displacement=True)

""" 
Generate a draft Kinetic Model
"""

# Generate the KineticModel

# Define the molecules that should be considered small-molecules
# These molecules will not be accounted explicitly in the kinetic mechanism as
# substrates and products
small_molecules = ['h_c', 'h_e']

model_gen = FromPyTFA(small_molecules=small_molecules)
this_skimpy_model = model_gen.import_model(this_pytfa_model,solution.raw)

export_to_yaml(this_skimpy_model, 'toy_draft.yml')
