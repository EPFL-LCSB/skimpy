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
from scipy.sparse import coo_matrix
from sympy import symbols
from sympy.utilities.autowrap import ufuncify


class ElasticityFunction:
    def __init__(self, expressions, variables,  parameters, shape):
        """
        Constructor for a precompiled function to compute elasticities
        numerically
        :param variables: a list of strings denoting
                                      the independent variables names
        :param expressions: dict of  non-zero sympy expressions for the rate of
                            change of a variable indexed by a tuple of the matrix position
                            e.g: (1,1)
        :param parameters:  list of parameter names
        :param shape: Tuple defining the over all matrix size e.g (10,30)

        """

        self.variables = variables
        self.expressions = expressions
        self.parameters = parameters
        self.shape = shape

        # Unpacking is needed as ufuncify only take ArrayTypes
        parameters = [x for x in self.parameters]
        variables = [x for x in variables]

        sym_vars = list(symbols(variables+parameters))

        # Awsome sympy magic
        # TODO problem with typs if any parameter ot variables is interpreted as interger
        # Make a function to compute every non zero entry in the matrix
        self.function = []
        self.coordinates = []
        for coord,exp in expressions.items():
            self.function.append(ufuncify(tuple(sym_vars),
                                         exp,
                                         backend='Cython'))
            self.coordinates.append(coord)


    def __call__(self, variables, parameters):
        """
        Return a sparse matrix type with elasticity values
        """
        parameter_values = [x for x in parameters.values()]

        input_vars = variables + parameter_values
        array_input = [array([input_var], dtype=double) for input_var in input_vars]

        values = [function(*array_input)[0] for function in self.function]
        rows, columns = zip(*self.coordinates)

        elasticiy_matrix = coo_matrix((values,
                                      (rows, columns)),
                                       shape=self.shape).tocsc()

        return elasticiy_matrix

