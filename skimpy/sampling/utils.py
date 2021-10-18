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
WITHOUT WARRANTIE CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""

import numpy as np
from sympy import Symbol

from numpy.linalg import eig as eigenvalues


def calc_max_eigenvalue(parameter_sample,
                        compiled_model,
                        concentration_dict,
                        flux_dict):

    """
    Sample one set of staturations using theano complied functions
    :param compiled_model:
    :param concentration_dict:
    :param flux_dict:
    :return:
    """
    reactions = compiled_model.reactions.values()
    fluxes = [flux_dict[this_reaction.name] for this_reaction in reactions]

    concentrations = np.array([concentration_dict[this_variable] for
                               this_variable in compiled_model.variables.keys()])

    # Check stability: real part of all eigenvalues of the jacobian is <= 0
    this_jacobian = compiled_model.jacobian_fun(fluxes, concentrations,
                                                parameter_sample)

    # largest_eigenvalue = eigenvalues(this_jacobian, k=1, which='LR',
    #                                 return_eigenvectors=False)
    # Test suggests that this is apparently much faster ....
    largest_eigenvalue = np.real(sorted(
        eigenvalues(this_jacobian.todense())[0]))[-1]

    return largest_eigenvalue


def calc_parameters( saturations,
                     compiled_model,
                     concentration_dict,
                     flux_dict,
                     parameters_to_resample=None,
                     fixed_parameters=None
                     ):

    symbolic_concentrations_dict = {Symbol(k):v
                                    for k,v in concentration_dict.items()}

    parameter_sample = {v.symbol: v.value for k,v in compiled_model.parameters.items()}

    # Update the concentrations which are parameters (Boundaries)
    parameters = compiled_model.parameters.keys()
    for k,v in symbolic_concentrations_dict.items():
        if str(k) in parameters:
            parameter_sample[k] = v


    #Set all vmax/flux parameters to 1.
    # TODO Generalize into Flux and Saturation parameters
    reactions = compiled_model.reactions.values()
    for this_reaction in reactions:
        vmax_param = this_reaction.parameters.vmax_forward
        parameter_sample[vmax_param.symbol] = 1.0

    if not hasattr(compiled_model,'saturation_parameter_function')\
       or not hasattr(compiled_model,'flux_parameter_function'):
        raise RuntimeError("Function for sampling not complied")


    # Calcualte the Km's
    compiled_model.saturation_parameter_function(
        saturations,
        parameter_sample,
        symbolic_concentrations_dict,
        parameters_to_resample,
        fixed_parameters
    )

    # Calculate the Vmax's
    compiled_model.flux_parameter_function(
        compiled_model,
        parameter_sample,
        symbolic_concentrations_dict,
        flux_dict
    )

    return parameter_sample