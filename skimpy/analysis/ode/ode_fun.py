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

from numpy import array, double
from sympy import symbols, Symbol

from skimpy.utils.compile_sympy import make_cython_function
from skimpy.utils.general import robust_index
from ...utils.tabdict import TabDict
from warnings import warn


class ODEFunction:
    def __init__(self, model, variables, expressions, parameters,
                 pool=None, with_time=False, custom_ode_update=None):
        """
        Constructor for a precompiled function to solve the ode epxressions
        numerically
        :param variables: a list of strings with variables names
        :param expressions: dict of sympy expressions for the rate of
                     change of a variable indexed by the variable name
        :param parameters: dict of parameters

        """
        self.variables = variables
        self.expressions = expressions
        self.model = model
        self.with_time = with_time
        self.custom_ode_update=custom_ode_update
        # Link to the model
        self._parameters = parameters

        # Unpacking is needed as ufuncify only take ArrayTypes
        the_param_keys = [x for x in self._parameters]
        the_variable_keys = [x for x in variables]

        if with_time:
            the_variable_keys = ['t',] + the_variable_keys

        sym_vars = list(symbols(the_variable_keys+the_param_keys))

        # Sort the expressions
        expressions = [self.expressions[x] for x in self.variables.values()]

        # Awsome magic
        self.function = make_cython_function(sym_vars, expressions, simplify=True, pool=pool)

    @property
    def parameters(self):
        model_params = self.model.parameters
        return TabDict((k, model_params[robust_index(k)].value)
                       for k in self._parameters)

    @parameters.setter
    def parameters(self, value):
        self._parameters = value

    def get_params(self):
        self._parameters_values = self.parameters.values()

    def __call__(self, t, y, ydot):
        if self.with_time:
            input_vars = [t,]+list(y)+list(self._parameters_values)
        else:
            input_vars = list(y)+list(self._parameters_values)
        self.function(input_vars, ydot)
        
        if not self.custom_ode_update is None:
            self.custom_ode_update( t, y, ydot)
