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
from numpy.random import sample

from sympy import symbols,Symbol
from skimpy.utils.compile_sympy import make_cython_function

class SaturationParameterFunction():
    def __init__(self,model,parameters,concentrations):


        sym_concentrations = [c for c in concentrations]

        self.sym_concentrations = sym_concentrations

        self.saturation_parameters = [v for k,v in parameters.items()
                                     if (v.hook is not None) \
                                     and (v.value is None) ]
        sym_saturations = []
        expressions = []
        for p in self.saturation_parameters:
            this_sat_symbol = Symbol("sigma_"+str(p.symbol))
            sym_saturations.append(this_sat_symbol )
            expressions.append((1-this_sat_symbol)*p.hook.symbol/this_sat_symbol)

        self.expressions = expressions
        self.sym_saturations = sym_saturations

        sym_vars = sym_saturations + sym_concentrations

        self.function = make_cython_function(sym_vars, expressions, simplify=False, pool=model.pool)

    def __call__(self,saturations,parameters,concentrations):
        _saturations = saturations
        _concentrations = np.array([concentrations[c] for c in self.sym_concentrations])

        input = np.concatenate((_saturations,_concentrations))
        saturation_parameter_values = np.zeros(len(self.saturation_parameters))

        self.function(input,saturation_parameter_values)

        # Assing saturation parameters
        for p,v in zip(self.saturation_parameters, saturation_parameter_values):
            parameters[p.symbol] = v
