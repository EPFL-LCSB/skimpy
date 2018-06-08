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
from numpy import array, double, reciprocal
from numpy import append as append_array
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import inv as sparse_inv
from sympy import symbols,Symbol
from sympy.printing.theanocode import theano_function

from skimpy.utils.tabdict import TabDict


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

        self.dummy = TabDict([])
        sym_vars = list(symbols(variables+parameters))

        # Awsome sympy magic
        # TODO problem with typs if any parameter ot variables is interpreted as interger
        # Make a function to compute every non zero entry in the matrix
        if 1 in expressions.values():
            dummy = Symbol('dummy_one')
            sym_vars += [dummy]
            self.dummy[dummy] = 1
            coordinates, expressions= zip(*[ (coord, expr*dummy) if expr == 1 else (coord, expr)
                                             for coord, expr in expressions.items()])
        else:
            coordinates, expressions = zip(*[ (coord, expr) for coord, expr in expressions.items()])

        rows, columns = zip(*coordinates)
        self.rows = rows
        self.columns = columns

        self.function = theano_function(sym_vars,expressions)


    def __call__(self, variables, parameters):
        """
        Return a sparse matrix type with elasticity values
        """
        parameter_values = array([parameters[x] for x in self.parameters.values()], dtype=double)

        dummy_values = array([x for x in self.dummy.values()])
        if not dummy_values:
            input_vars = append_array(variables, parameter_values)
        else:
            input_vars = append_array(variables , parameter_values ,dummy_values)

        values = self.function(*input_vars)

        elasticiy_matrix = coo_matrix((values,
                                      (self.rows, self.columns)),
                                       shape=self.shape).tocsc()

        return elasticiy_matrix

    def get_dependent_weights(self, concentration_vector, L0, all_independent_ix, all_dependent_ix):

        # The dependent weights have dimensions of moieties x independent metabolites
        # The current L0 gives the relation L0*[xi|xd] = C

        # Concentrations
        X = concentration_vector
        Xi = X[all_independent_ix]
        Xd = X[all_dependent_ix]

        # L0 = [Fd | Fi]
        # Thus Fd*x_d = -Fi*xi + C
        # xd = (Fd^-1).(-Fi*xi + C)
        # dxd/dxi = (Fd^-1).(-Fi)
        # The dependent weights are Qd :
        # Qd = dln(xd)/dln(xi)
        # Qd = dxd / xd * xi / dxi
        # Qd = ( xd^-1 ) * dxd/dxi * xi

        # Fi Factors for in dependent concentrations
        Fi = L0[:, all_independent_ix]
        # Fd Factors for dependent concentrations
        Fd = L0[:, all_dependent_ix]
        # Qd the scaling matrix
        dxd_dxi = sparse_inv(Fd).dot(-1*Fi)
        Qd = reciprocal(Xd).dot(dxd_dxi).dot(Xi)

        return Qd

