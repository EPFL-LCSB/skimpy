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
from sympy import symbols, Array
from sympy.utilities.autowrap import ufuncify


class ODEFunction:
    def __init__(self, variables, expr, parameters):
        self.variables  = variables
        self.expr = expr
        self.parameters = parameters

        # Unpacking is needed as ufuncify only take ArrayTypes
        the_param_keys = [x for x in self.parameters]
        the_variable_keys = [x for x in variables]
        sym_vars = list(symbols(the_variable_keys+the_param_keys))

        # Sort the expressions
        the_expressions = [self.expr[x] for x in self.variables]

        # Awsome sympy magic
        # TODO problem with typs if any parameter ot vairabls is interpreted as interger
        self.function = []
        for exp in the_expressions:
           self.function.append(ufuncify(tuple(sym_vars),
                                         exp,
                                         backend='Cython'))

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self,value):
        """
        Would-be optimization hack to avoid looking up thr whole dict at each
        iteration step in __call__

        :param value:
        :return:
        """
        self._parameters = value
        self.parameter_values = [x for x in self.parameters.values()]


    def __call__(self,t,y):
        input_vars = list(y)+self.parameter_values
        #result = self.functin
        array_input = [np.array([input_var])  for input_var in  input_vars  ]
        results = [f(*array_input) for f in self.function ]
        return np.array(results)
