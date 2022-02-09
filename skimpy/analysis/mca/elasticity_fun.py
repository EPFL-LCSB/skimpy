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
from numpy import array, double, reciprocal,zeros
from numpy import append as append_array

# Test wise
from scipy.sparse import SparseEfficiencyWarning
import warnings
warnings.simplefilter('ignore',SparseEfficiencyWarning)

from scipy.sparse import coo_matrix
from scipy.sparse import diags, find
from scipy.sparse.linalg import inv as sparse_inv
from sympy import symbols,Symbol

from skimpy.utils.tabdict import TabDict
from skimpy.utils.compile_sympy import make_cython_function
from skimpy.utils.general import robust_index

class ElasticityFunction:
    def __init__(self, expressions, respective_variables, variables,  parameters, shape, pool=None):
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
        self.respective_variables = respective_variables
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

        coordinates, expressions= zip(*[ (coord, expr) for coord, expr in expressions.items()])

        rows, columns = zip(*coordinates)
        self.rows = rows
        self.columns = columns

        # self.function = theano_function(sym_vars, expressions,
        #                                 on_unused_input='ignore')

        self.function = make_cython_function(sym_vars, expressions, pool=pool,
                                             simplify=True, optimize=True)

    def __call__(self, variables, parameters):
        """
        Return a sparse matrix type with elasticity values
        """
        parameter_values = array([parameters[x] for x in
                                  self.parameters.values()], dtype=double)

        input_vars = append_array(variables , parameter_values)

        values = array(zeros(len(self.expressions)),dtype=double)

        self.function(input_vars, values)

        elasticiy_matrix = coo_matrix((values,
                                      (self.rows, self.columns)),
                                       shape=self.shape).tocsc()

        return elasticiy_matrix

    def get_dependent_weights(self, concentration_vector,
                              L0,
                              all_independent_ix,
                              all_dependent_ix,
                              volume_ratios = None):

        # TODO This derivation does not allow cross dependencies of dependent metabolites!
        # Usually you can find basis that omit this

        # The dependent weights have dimensions of moieties x independent metabolites
        # The current L0 gives the relation L0*[xi|xd] = C

        # Concentrations
        X = array(concentration_vector)
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

        if volume_ratios is None:
            # Fi Factors for in dependent concentrations
            Fi = L0[:, all_independent_ix]
            # Fd Factors for dependent concentrations
            Fd = L0[:, all_dependent_ix]

        else:
            v_d_ = diags( reciprocal(volume_ratios[all_dependent_ix])).tocsc()
            Fd = L0[:, all_dependent_ix].dot(v_d_)
            v_i_ = diags( reciprocal(volume_ratios[all_independent_ix])).tocsc()
            Fi = L0[:, all_independent_ix].dot(v_i_)

        dxd_dxi = sparse_inv(Fd).dot(-Fi)
        if len(dxd_dxi.shape) == 1:
            dxd_dxi = dxd_dxi[0]
        # Qd = dxd_dxi.multiply(Xi).T.multiply(reciprocal(Xd)).T
        XD = diags(reciprocal(Xd), 0).tocsc()
        XI = diags(Xi, 0).tocsc()
        Qd = XD.dot(dxd_dxi).dot(XI)

        return Qd

