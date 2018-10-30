
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

def add_min_log_displacement(tmodel,min_log_displacement, inplace=True):
    if inplace:
        temp_model = tmodel
    else:
        temp_model = tmodel.copy()
        temp_model.repair()

    for ln_gamma in temp_model.thermo_displacement:
         if ln_gamma.variable.primal > 0:
             ln_gamma.variable.lb = min_log_displacement
         if ln_gamma.variable.primal < 0:
             ln_gamma.variable.ub = -min_log_displacement

    return temp_model
