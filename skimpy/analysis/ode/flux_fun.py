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

import numpy as np
from sympy import symbols
from skimpy.utils.compile_sympy import make_cython_function


class FluxFunction:
    def __init__(self, variables, expr, parameters, pool=None):
        """
        Constructor for a precompiled function to solve the ode epxressions
        numerically
        :param variables: a list of strings with variables names
        :param expr: dict of sympy expressions for the rate of
                     change of a variable indexed by the variable name
        :param parameters: dict of parameters with parameter values

        """
        self.variables = variables
        self.expr = expr
        self.parameters = parameters

        # Unpacking is needed as ufuncify only take ArrayTypes
        the_param_keys = [x for x in self.parameters]
        the_variable_keys = [x for x in variables]
        sym_vars = list(symbols(the_variable_keys+the_param_keys))

        self.function = make_cython_function(sym_vars, expr.values(), simplify=True, pool=pool)


    def __call__(self,concentrations,  parameters=None):
        # Todo handle different input types
        variables = [concentrations[str(x)] for x in self.variables]

        if parameters is None:
            input_vars = list(variables)+list(self.parameters.values())
        else:
            input_vars = list(variables) \
                         + [parameters[x] for x in self.parameters]

        fluxes = np.zeros(len(self.expr))

        self.function(input_vars, fluxes)

        return {k:v for k,v in zip(list(self.expr.keys()) , fluxes)}
