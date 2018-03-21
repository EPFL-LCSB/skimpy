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

from abc import ABC, abstractmethod
from collections import namedtuple
from numpy import real
from numpy import all as np_all
from numpy.random import sample
from scipy.sparse.linalg import eigs as eigenvalues
from sympy import sympify

class ParameterSampler(ABC):
    def __init__(self, parameters=None):
        """

        :param parameters:
        """
        self.parameters = parameters

    @property
    @abstractmethod
    def Parameters(self):
        """
        Parameter type specified for the parmeters samples
        :return:
        """


    @abstractmethod
    def sample(self):
        """

        :return:
        """


class SimpleParameterSampler(ParameterSampler):
    """
    A simple parameter sampler that samples stable model parameters
    with respect to a steady state flux and concentration state
    """

    Parameters = namedtuple('Parameters', ['n_samples'])
    # TODO Talk to Pierre / Misko about simple sampler parameters
    # if parameters are not defined put default values
    Parameters.__new__.__defaults__ = (None,) * len(Parameters._fields)


    def sample(self,compiled_model,flux_dict,concentration_dict):

        parameter_population = []

        # Unpack fluxes and concentration into arrays consitent with the
        # compiled functions

        fluxes = [flux_dict[this_reaction.name] for this_reaction in compiled_model.reactions]
        concentrations = [concentration_dict[this_variable] for this_variable in compiled_model.variables.keys()]

        for this_sample in range(self.parameters.n_samples):
            parameter_sample = {}
            # Sample parameters for every reaction
            for this_reaction in compiled_model.reactions:
                this_parameters = {'k_equilibrium_'+this_reaction.name: this_reaction.mechanism.parameters.k_equilibrium,
                                   'vmax_forward_'+this_reaction.name: 1.0}

                # Loop over the named tuple
                for this_name, this_reactant in this_reaction.mechanism.substrates._asdict().iteritems():
                    # Sample a saturation
                    this_saturation = sample()
                    this_concentration = concentration_dict[this_reactant]
                    this_km_name  = 'km_'+'_'+this_name+'_'+this_reactant.name+'_'+this_reaction.name
                    this_parameters[this_km_name] = (1.0-this_saturation) * this_concentration / this_saturation


                # Calculate the vmax
                this_net_reaction_rate = this_reaction.mechanism.reaction_rates['v_net']
                this_parameter_subs = [(sympify(var), val) for var, val in this_parameters.items()]
                normed_net_reaction_rate = this_net_reaction_rate.subs(this_parameter_subs).evalf()
                this_vmax = flux_dict[this_reaction.name]/normed_net_reaction_rate

                this_parameters['vmax_forward_'+this_reaction.name] = this_vmax

                # Update the dict with explicit model parameters
                parameter_sample.update(this_parameters)
                # TODO would thi be nicer???
                # parameter_sample[this_reaction.name] = this_reaction.mechanism.Parameters(**this_parameters)

            concentrations = [concentration_dict[this_variable] for this_variable in compiled_model.variables.keys()]

            # Parametrize the Jacobian function
            compiled_model.jacobian_fun.parameters = parameter_sample

            # Check stability: real part of all eigenvalues of the jacobian is <= 0
            this_jacobian = compiled_model.jacobian_fun(fluxes, concentrations)
            largest_eigenvalue = eigenvalues(this_jacobian,k = 1, wich ='LR',  return_eigenvectors=False )
            is_stable = largest_eigenvalue < 0

            if is_stable:
                parameter_population.append(parameter_sample)

        return parameter_population






