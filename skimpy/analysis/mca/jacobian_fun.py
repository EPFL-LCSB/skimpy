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

from numpy import array

from scipy.sparse import diags
from scipy.sparse.linalg import inv as sparse_inv


class JacobianFunction:
    def __init__(self,
                 reduced_stoichometry,
                 independent_elasticity_function,
                 depednent_elasticity_function,
                 weights_dependent_metabolites
                 ):

        self.reduced_stoichometry = reduced_stoichometry
        self.depednent_elasticity_function  = depednent_elasticity_function
        self.independent_elasticity_function = independent_elasticity_function
        self.weights_dependent_metabolites = weights_dependent_metabolites

    def __call__(self, fluxes, concentrations, parameters):

        #Calculate the Jacobian
        flux_matrix = diags(array(fluxes), 0)
        concentration_matrix = diags(array(concentrations), 0)

        inv_concentration_matrix = sparse_inv(concentration_matrix)

        # Elasticity matrix
        if len(self.weights_dependent_metabolites) == 0:
            elasticity_matrix = self.independent_elasticity_function(concentrations,parameters)
        else:
            elasticity_matrix = self.independent_elasticity_function(concentrations, parameters)
            elasticity_matrix += self.depednent_elasticity_function(concentrations, parameters)\
                                 .dot(self.weights_dependent_metabolites)

        jacobian = self.reduced_stoichometry.dot(flux_matrix)\
                    .dot(elasticity_matrix)\
                    .dot(inv_concentration_matrix)

        return jacobian
