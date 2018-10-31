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

from numpy import array, zeros

from scipy.sparse import diags
from scipy.sparse.linalg import inv as sparse_inv


class ConcentrationControlFunction:
    def __init__(self,
                 reduced_stoichometry,
                 independent_elasticity_function,
                 dependent_elasticity_function,
                 parameter_elasticity_function,
                 conservation_relation,
                 independent_variable_ix,
                 dependent_variable_ix,
                 ):

        self.reduced_stoichometry = reduced_stoichometry
        self.dependent_elasticity_function = dependent_elasticity_function
        self.independent_elasticity_function = independent_elasticity_function
        self.parameter_elasticity_function = parameter_elasticity_function
        self.independent_variable_ix = independent_variable_ix
        self.dependent_variable_ix = dependent_variable_ix
        self.conservation_relation = conservation_relation

    def __call__(self, fluxes, concentrations, parameter_population):

        # Calculate the Concentration Control coefficients
        # Log response of the concentration with respect to the log change in a Parameter
        #
        # C_Xi_P = -(N_r*V*E_i + N_r*V*E_d*Q_i)(N_r*V*Pi)
        #
        num_parameters = len(self.parameter_elasticity_function.expressions)
        num_concentration = len(concentrations)
        population_size = len(parameter_population)

        concentration_control_coefficients = zeros((num_concentration,num_parameters,population_size))

        for i,parameters in enumerate(parameter_population):

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

            N_E_V = self.reduced_stoichometry.dot(flux_matrix).dot(elasticity_matrix)
            N_E_V_inv = sparse_inv(N_E_V)

            parameter_elasticity_matrix = self.parameter_elasticity_function(concentrations, parameters)

            N_E_P = self.reduced_stoichometry.dot(flux_matrix).dot(parameter_elasticity_matrix)

            this_cc = - N_E_V_inv.dot(N_E_P)
            concentration_control_coefficients[:,:,i] = this_cc.todense()


        return concentration_control_coefficients
