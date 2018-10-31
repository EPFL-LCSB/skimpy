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


class FluxControlFunction:
    def __init__(self,
                 reduced_stoichometry,
                 independent_elasticity_function,
                 dependent_elasticity_function,
                 parameter_elasticity_function,
                 conservation_relation,
                 independent_variable_ix,
                 dependent_variable_ix,
                 concentration_control_fun,
                 ):

        self.reduced_stoichometry = reduced_stoichometry
        self.dependent_elasticity_function = dependent_elasticity_function
        self.independent_elasticity_function = independent_elasticity_function
        self.parameter_elasticity_function = parameter_elasticity_function
        self.independent_variable_ix = independent_variable_ix
        self.dependent_variable_ix = dependent_variable_ix
        self.conservation_relation = conservation_relation

        self.concentration_control_fun = concentration_control_fun

    def __call__(self, fluxes, concentrations, parameters):

        # Calculate the Flux Control coefficients
        # Log response of the concentration with respect to the log change in a Parameter
        #
        # C_V_P = (E_i + E_d*Q_i)*C_Xi_P + Pi
        #

        flux_matrix = diags(array(fluxes), 0).tocsc()

        # Elasticity matrix

        if self.conservation_relation.nnz == 0:
            # If there are no moieties
            elasticity_matrix = self.independent_elasticity_function(concentrations,parameters)

        else:
            # If there are moieties
            ix = self.independent_variable_ix
            elasticity_matrix = self.independent_elasticity_function(concentrations, parameters)
            dependent_weights = self.dependent_elasticity_function.\
                get_dependent_weights(
                                concentration_vector=concentrations,
                                L0=self.conservation_relation,
                                all_dependent_ix=self.dependent_variable_ix,
                                all_independent_ix=self.independent_variable_ix,
                            )

            # Calculate the effective elasticises
            elasticity_matrix += self.dependent_elasticity_function(concentrations, parameters)\
                                 .dot(dependent_weights)

        C_Xi_P = self.concentration_control_fun(fluxes, concentrations, parameters)

        parameter_elasticity_matrix = self.parameter_elasticity_function(concentrations, parameters)

        flux_control_coefficients = elasticity_matrix.dot(C_Xi_P) + parameter_elasticity_matrix

        return flux_control_coefficients
