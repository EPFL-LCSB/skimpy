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
from ...utils.tabdict import TabDict
from warnings import warn

def robust_index(in_var):
    """
    Indexing can be done with symbols or strings representing the symbol,
    so we harmonize it by returning the name of the symbol if the input is of
    type symbol

    :param in_var:
    :type in_var: str or sympy.Symbol
    :return:
    """

    if isinstance(in_var, str):
        return in_var
    elif isinstance(in_var, Symbol):
        return in_var.name
    else:
        raise TypeError('Value should be of type str or sympy.Symbol')

class ODEFunction:
    def __init__(self, model, variables, expressions, parameters, pool=None):
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
        # self._parameter_values = TabDict([])

        # Link to the model
        self._parameters = parameters

        # Set values
        for k,v in parameters.items():
            self.parameters[robust_index(k)] = v

        # Unpacking is needed as ufuncify only take ArrayTypes
        the_param_keys = [x for x in self._parameters]
        the_variable_keys = [x for x in variables]
        sym_vars = list(symbols(the_variable_keys+the_param_keys))

        # Sort the expressions
        expressions = [self.expressions[x] for x in self.variables.values()]

        # Awsome magic
        self.function = make_cython_function(sym_vars, expressions, simplify=False, pool=pool)

    @property
    def parameters(self):
        return TabDict((k,self.model.parameters[robust_index(k)].value)
                                   for k in self._parameters)

    @parameters.setter
    def parameters(self, value):
        self._parameters = value

    # @property
    # def parameter_values(self):
    #     # if not self._parameter_values:
    #     #     raise Exception('No parameters have been set')
    #     # else:
    #     return TabDict((k,self.model.parameters[robust_index(k)].value)
    #                    for k in self.parameters)
    #
    # @parameter_values.setter
    # def parameter_values(self,value):
    #     """
    #     Would-be optimization hack to avoid looking up thr whole dict at each
    #     iteration step in __call__
    #
    #     :param value:
    #     :return:
    #     """
    #     #self._parameters = value
    #     # self._parameter_values = [value[x] for x in self.parameters.values()]
    #
    #     for k,v in value.items():
    #         if v is None:
    #             # No assignment is to be done
    #             continue
    #
    #         try:
    #             self.parameters[robust_index(k)].value = v
    #         except KeyError:
    #             # raise KeyError('Parameter is not in the model.')
    #             warn('Tried to assign a value to parameter {}. '
    #                  'Parameter is not in the model'.format(k))
    #             self.parameters[robust_index(k)] = v

    def __call__(self, t, y, ydot):
        input_vars = list(y)+list(self.parameters.values())
        self.function(input_vars,ydot)
