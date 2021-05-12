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
WITHOUT WARRANTIE CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""
import numpy as np

from sympy import symbols,Symbol
from skimpy.utils.compile_sympy import make_cython_function

class FluxParameterFunction():
    def __init__(self,
                 model,
                 parameters,
                 concentration_dict):

        self.sym_concentrations = [c for c in concentration_dict]

        self.sym_parameters = [p.symbol for p in parameters.values()
                               if p.symbol not in self.sym_concentrations ]

        self.expressions = [ rxn.mechanism.reaction_rates['v_net']
                             for rxn in model.reactions.values()]

        sym_vars = self.sym_parameters+self.sym_concentrations
        self.function = make_cython_function(sym_vars, self.expressions, simplify=True,
                                             pool=model.pool)

    def __call__(self,
                 model,
                 parameters,
                 concentration_dict,
                 flux_dict):
        _parameters = [parameters[p] for p in self.sym_parameters]
        _concentrations = [concentration_dict[c] for c in self.sym_concentrations]

        input = _parameters + _concentrations
        flux_parameter_values = np.zeros(len(model.reactions))

        self.function(input,flux_parameter_values)

        _fluxes = np.array([flux_dict[rxn.name] for rxn in model.reactions.values() ])
        flux_parameter_values = _fluxes / flux_parameter_values

        if np.any(flux_parameter_values < 0):
            ixs = np.where(flux_parameter_values < 0)[0]
            model.logger.info("Fluxes {} are not aligned with deltaG values!".format([model.reactions.iloc(i)[0] for i in ixs]))
            raise ValueError

        for rxn,v in zip(model.reactions.values(),flux_parameter_values):
            parameters[rxn.parameters.vmax_forward.symbol] = v
