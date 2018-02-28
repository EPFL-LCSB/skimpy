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

from sympy import symbols, Array
from sympy.utilities.autowrap import ufuncify



class ElasticityFunction:
    def _init_(self,
               independent_variables,
               dependent_variables,
               expressions,
               parameters):
        """
        Constructor for a precompiled function to compute elasticities
        numerically
        :param independent_variables: a list of strings denoting
                                      the independent variables names
        :param dependent_variables: a list of strings denoting
                                      the dependent variables names
        :param expressions: dict of sympy expressions for the rate of
                            change of a variable indexed by the variable
                            name
        :param parameters: dict of parameters with parameter values

        """
        self.independent_variables = independent_variables
        self.dependent_variables = dependent_variables
        self.expressions = expressions
        self.parameters = parameters


        # Unpacking is needed as ufuncify only take ArrayTypes
        the_param_keys = [x for x in self.parameters]
        the_independent_variable_keys = [x for x in independent_variables]
        the_dependent_variable_keys = [x for x in dependent_variables]

        sym_vars = list(symbols(independent_variables+dependent_variables+the_param_keys))

        # Sort the expressions
        the_expressions = [self.expressions[x] for x in self.the_independent_variable_keys]

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

    def __call__(self,
                 independent_variables,
                 dependent_variables,
                 expressions,
                 parameters):
        pass
