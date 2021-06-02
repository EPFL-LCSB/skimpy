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
from numpy import array, zeros

from scipy.sparse import diags
from scipy.sparse.linalg import inv as sparse_inv

from skimpy.utils.tensor import Tensor
from skimpy.utils.namespace import SPLIT, NET


class FluxControlFunction:
    def __init__(self,
                 model,
                 reduced_stoichometry,
                 independent_elasticity_function,
                 dependent_elasticity_function,
                 parameter_elasticity_function,
                 conservation_relation,
                 independent_variable_ix,
                 dependent_variable_ix,
                 concentration_control_fun,
                 mca_type=NET,
                 ):

        self.model = model
        self.reduced_stoichometry = reduced_stoichometry
        self.dependent_elasticity_function = dependent_elasticity_function
        self.independent_elasticity_function = independent_elasticity_function
        self.parameter_elasticity_function = parameter_elasticity_function
        self.independent_variable_ix = independent_variable_ix
        self.dependent_variable_ix = dependent_variable_ix
        self.conservation_relation = conservation_relation

        self.concentration_control_fun = concentration_control_fun

        self.mca_type = mca_type


    def __call__(self, flux_dict, concentration_dict, parameter_population):

        # Calculate the Flux Control coefficients
        # Log response of the concentration with respect to the log change in a Parameter
        #
        # C_V_P = (E_i + E_d*Q_i)*C_Xi_P + Pi
        #

        fluxes = [flux_dict[r] for r in self.model.reactions]
        concentrations = [concentration_dict[r] for r in self.model.reactants]

        num_parameters = len(self.parameter_elasticity_function.respective_variables)

        if self.mca_type == NET:
            num_fluxes = len(fluxes)
        elif self.mca_type == SPLIT:
            num_fluxes = len(fluxes)*2

        population_size = len(parameter_population)


        flux_control_coefficients = zeros((num_fluxes,num_parameters,population_size))

        for i, parameters in enumerate(parameter_population):

            if self.concentration_control_fun.volume_ratio_function is None:
                volume_ratios = array([1, ] * len(concentrations) )
            else:
                volume_ratios = self.concentration_control_fun.volume_ratio_function(parameters)

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
                                    volume_ratios=volume_ratios
                                )

                # Calculate the effective elasticises
                elasticity_matrix += self.dependent_elasticity_function(concentrations, parameters)\
                                     .dot(dependent_weights)

            C_Xi_P = self.concentration_control_fun(flux_dict, concentration_dict, [parameters])._data

            parameter_elasticity_matrix = self.parameter_elasticity_function(concentrations, parameters)

            this_cc = elasticity_matrix.dot(C_Xi_P[:,:,0]) + parameter_elasticity_matrix
            flux_control_coefficients[:,:,i] = this_cc

        if self.mca_type == NET:
            flux_index = pd.Index(self.model.reactions.keys(), name="flux")
        elif self.mca_type == SPLIT:
            fwd_fluxes = ["fwd_"+k for k in self.model.reactions.keys()]
            bwd_fluxes = ["bwd_"+k for k in self.model.reactions.keys()]
            flux_index = pd.Index(fwd_fluxes+bwd_fluxes, name="flux")

        parameter_index = pd.Index(self.parameter_elasticity_function.respective_variables, name="parameter")
        sample_index = pd.Index(range(population_size), name="sample")

        tensor_fcc = Tensor(flux_control_coefficients, [flux_index,parameter_index,sample_index])

        return tensor_fcc
