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
        for coord, exp in expressions.items():
            this_sym_vars = exp.free_symbols
            this_sym_var_ix = [i for i,e in enumerate(sym_vars) if e in this_sym_vars]
            this_ordered_sym_vars = [e for i, e in enumerate(sym_vars) if e in this_sym_vars]
            # Cast a dummy value for elasticities = 1
            if not this_ordered_sym_vars:
                this_ordered_sym_vars = sym_vars[0:2]
                this_sym_var_ix = [0,1]

            self.function.append((ufuncify(tuple(this_ordered_sym_vars),
                                         exp,
                                         backend='Cython'),this_sym_var_ix))
            self.coordinates.append(coord)


    def __call__(self, variables, parameters):
        """
        Return a sparse matrix type with elasticity values
        """
        parameter_values = array([parameters[x] for x in self.parameters.values()], dtype=double)

        input_vars = append_array(variables , parameter_values)
        array_input = array([array([input_var], dtype=double) for input_var in input_vars])
        values = [function(*array_input[ix])[0] for function, ix in self.function]
        rows, columns = zip(*self.coordinates)

        elasticiy_matrix = coo_matrix((values,
                                      (rows, columns)),
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

