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

from numpy import array, zeros, double
from numpy import append as append_array
from sympy import symbols
from sympy import diff

from scipy.sparse import coo_matrix

from skimpy.utils.compile_sympy import make_cython_function
from skimpy.utils.general import join_dicts


class SymbolicJacobianFunction:

    def __init__(self, variables, ode_expressions, parameters, pool=None):
        """
        Constructor for a precompiled function to compute epxressions
        numerically
        :param variables: a list of strings with variables names
        :param expr: dict of sympy expressions for the rate of
                     change of a variable indexed by the variable name
        :param parameters: dict of parameters

        """
        self.variables = variables
        self.ode_expressions = ode_expressions
        self.parameters = parameters
        self.shape = (len(variables), len(ode_expressions) )

        # Unpacking is needed as ufuncify only take ArrayTypes
        parameters = [x for x in self.parameters]
        variables = [x for x in variables]

        sym_vars = list(symbols(variables+parameters))

        # Awsome sympy magic
        # Awsome sympy magi
        # TODO problem with typs if any parameter ot variables is interpreted as interger
        # Make a function to compute every non zero entry in the matrix

        # Compute the Jacobian
        self.expressions = make_symbolic_jacobian(self.variables.values(),ode_expressions, pool=pool)


        coordinates, expressions= zip(*[ (coord, expr) for coord, expr in self.expressions.items()])

        rows, columns = zip(*coordinates)
        self.rows = rows
        self.columns = columns

        # self.function = theano_function(sym_vars, expressions,
        #                                 on_unused_input='ignore')

        self.function = make_cython_function(sym_vars, expressions, pool=pool, simplify=False)

    def __call__(self, fluxes, concentrations, parameters):
        """
        Return a sparse matrix type with elasticity values
        """
        parameter_values = array([parameters[x.symbol] for x in self.parameters.values()], dtype=double)

        input_vars = append_array(concentrations , parameter_values)

        values = array(zeros(len(self.expressions)),dtype=double)

        self.function(input_vars, values)

        jacobian = coo_matrix((values,
                              (self.rows, self.columns)),
                               shape=self.shape).tocsc()

        return jacobian


def make_symbolic_jacobian(variables,ode_expressions, pool=None):
    # List of Vars and ode_

    if pool is None:
        expressions = {}
        for i, var_i in enumerate(variables):
            for j, var_j in enumerate(variables):
                derivative = diff(ode_expressions[var_j], var_i)
                if derivative != 0:
                    expressions[(i,j)] = derivative
    else:

        pickled_variables = [v for v in variables]
        pickled_ode_expressions = {k:v for k,v in ode_expressions.items()}

        inputs = [(i,var_i,pickled_variables,pickled_ode_expressions)
                  for i,var_i in enumerate(variables) ]

        row_slices = pool.map(make_symbolic_jacobian_row, inputs)

        expressions = join_dicts(row_slices)

    return expressions



def make_symbolic_jacobian_row(input):
    i, var_i, variables, ode_expressions = input
    expressions_row_slice = {}
    for j, var_j in enumerate(variables):
        derivative = diff(ode_expressions[var_j], var_i)
        if derivative != 0:
            expressions_row_slice[(i, j)] = derivative

    return expressions_row_slice
