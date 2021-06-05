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

import pandas as pd
from numpy import array, zeros, append

from scipy.sparse import diags
from scipy.sparse.linalg import inv as sparse_inv
from scipy.sparse import hstack

from skimpy.utils.tensor import Tensor
from skimpy.utils.namespace import SPLIT, NET
from skimpy.analysis.mca.utils import get_reversible_fluxes

class ConcentrationControlFunction:
    def __init__(self,
                 model,
                 reduced_stoichometry,
                 independent_elasticity_function,
                 dependent_elasticity_function,
                 parameter_elasticity_function,
                 volume_ratio_function,
                 conservation_relation,
                 independent_variable_ix,
                 dependent_variable_ix,
                 mca_type=NET,
                 displacement_function=None,
                 ):

        self.model = model
        self.reduced_stoichometry = reduced_stoichometry
        self.dependent_elasticity_function = dependent_elasticity_function
        self.independent_elasticity_function = independent_elasticity_function
        self.parameter_elasticity_function = parameter_elasticity_function
        self.independent_variable_ix = independent_variable_ix
        self.dependent_variable_ix = dependent_variable_ix
        self.conservation_relation = conservation_relation
        self.volume_ratio_function = volume_ratio_function

        self.displacement_function = displacement_function
        self.mca_type = mca_type

    def __call__(self,  flux_dict, concentration_dict, parameter_population):

        # Calculate the Concentration Control coefficients
        # Log response of the concentration with respect to the log change in a Parameter
        #
        # C_Xi_P = -(N_r*V*E_i + N_r*V*E_d*Q_i)(N_r*V*Pi)
        #

        fluxes = [flux_dict[r] for r in self.model.reactions]

        # Only consider independent concentrations
        concentrations = [concentration_dict[r] for r in self.model.reactants]

        num_parameters = len(self.parameter_elasticity_function.respective_variables)
        num_concentration = len(self.independent_variable_ix)
        population_size = len(parameter_population)

        concentration_control_coefficients = zeros((num_concentration, num_parameters, population_size))


        for i,parameters in enumerate(parameter_population):

            # net control coeffecients
            if self.mca_type == NET:
                flux_matrix = diags(array(fluxes), 0).tocsc()
                effective_reduced_stoichiometry = self.reduced_stoichometry

            # fwd and bwd control coeffecients
            elif self.mca_type == SPLIT:
                displacements = self.displacement_function(concentration_dict, parameters=parameters)

                forward_fluxes, backward_fluxes = get_reversible_fluxes(flux_dict,
                                                                        displacements,
                                                                        self.model.reactions)

                flux_matrix = diags(append(forward_fluxes, backward_fluxes), 0).tocsc()
                effective_reduced_stoichiometry = hstack([self.reduced_stoichometry, -self.reduced_stoichometry],
                                                         format='csc')

            if self.volume_ratio_function is None:
                volume_ratios = array([1, ] * len(concentrations) )
            else:
                volume_ratios = self.volume_ratio_function(parameters)

            # Elasticity matrix

            if self.conservation_relation.nnz == 0:
                # If there are no moieties
                volume_ratio_matrix = diags(array(volume_ratios)).tocsc()

                elasticity_matrix = self.independent_elasticity_function(concentrations,parameters)

            else:
                # If there are moieties
                ix = self.independent_variable_ix

                volume_ratio_matrix = diags(array(volume_ratios)[ix]).tocsc()

                elasticity_matrix = self.independent_elasticity_function(concentrations, parameters)

                dependent_weights = self.dependent_elasticity_function.\
                    get_dependent_weights(
                                    concentration_vector=concentrations,
                                    L0=self.conservation_relation,
                                    all_dependent_ix=self.dependent_variable_ix,
                                    all_independent_ix=self.independent_variable_ix,
                                    volume_ratios=volume_ratios
                                )

                # Calculate the effective elasticises
                elasticity_matrix += self.dependent_elasticity_function(concentrations, parameters)\
                                     .dot(dependent_weights)



            N_E_V = volume_ratio_matrix.dot(effective_reduced_stoichiometry).dot(flux_matrix).dot(elasticity_matrix)
            N_E_V_inv = sparse_inv(N_E_V)

            parameter_elasticity_matrix = self.parameter_elasticity_function(concentrations, parameters)

            N_E_P = volume_ratio_matrix.dot(effective_reduced_stoichiometry).dot(flux_matrix).dot(parameter_elasticity_matrix)

            this_cc = - N_E_V_inv.dot(N_E_P)
            concentration_control_coefficients[:,:,i] = this_cc.todense()

        concentration_index = pd.Index([self.model.reactants.iloc(i)[0] for i in self.independent_variable_ix],
                                       name="concentration")
        parameter_index = pd.Index(self.parameter_elasticity_function.respective_variables, name="parameter")
        sample_index = pd.Index(range(population_size), name="sample")

        tensor_ccc = Tensor(concentration_control_coefficients, [concentration_index,parameter_index,sample_index])

        return tensor_ccc

