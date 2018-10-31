# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2018 Laboratory of Computational Systems Biotechnology (LCSB),
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

from pytfa.analysis import  variability_analysis,           \
                            apply_reaction_variability,     \
                            apply_generic_variability

from pytfa.optim.variables import DeltaG,DeltaGstd,LogConcentration,ThermoDisplacement
from pytfa.analysis import sample

thermo_vars = [DeltaGstd,LogConcentration,ThermoDisplacement]


def get_reference(tmodel, num_samples, n_processes=1 ):
    pass


def _sample(tmodel, num_samples, n_processes=1):

    tva_fluxes = variability_analysis(tmodel, kind='reactions')
    tva_thermo = variability_analysis(tmodel, kind=thermo_vars)

    tight_model = apply_reaction_variability(tmodel, tva_fluxes)
    tight_model = apply_generic_variability(tight_model, tva_thermo)

    sampling = sample(tight_model, num_samples, processes=n_processes)

    return sampling




