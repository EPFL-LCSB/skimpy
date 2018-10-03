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

from sympy import symbols,Symbol
from sympy.printing.theanocode import theano_function

class FluxParameterFunction():
    def __init__(self,
                 model,
                 parameters,
                 concentration_dict,
                 flux_dict):

        self.sym_concentrations = [c for c in concentration_dict]
        self.sym_fluxes = [Symbol(c) for c in flux_dict]

        self.sym_parameters = [p.symbol for p in parameters.values()
                               if p.symbol not in self.sym_concentrations ]

        self.expressions = [ Symbol(rxn.name) / rxn.mechanism.reaction_rates['v_net']
                             for rxn in model.reactions.values()]

        sym_vars = self.sym_parameters+self.sym_concentrations+self.sym_fluxes
        self.function = theano_function(sym_vars, self.expressions,
                                        on_unused_input='ignore')

    def __call__(self,
                 model,
                 parameters,
                 concentration_dict,
                 flux_dict):
        _parameters = [parameters[p] for p in self.sym_parameters]
        _concentrations = [concentration_dict[c] for c in self.sym_concentrations]
        _fluxes = [flux_dict[str(c)] for c in self.sym_fluxes]

        input = _parameters + _concentrations + _fluxes
        flux_parameter_values = self.function(*input)

        for rxn,v in zip(model.reactions.values(),flux_parameter_values):
            parameters[rxn.parameters.vmax_forward.symbol] = v
