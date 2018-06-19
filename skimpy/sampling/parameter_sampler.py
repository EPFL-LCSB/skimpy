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

import numpy as np
from numpy.random import sample
from scipy.sparse.linalg import eigs as eigenvalues
from sympy import sympify, Symbol

from skimpy.utils.namespace import *

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
        Parameter type specified for the parameters samples
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

    def sample(self, compiled_model, flux_dict, concentration_dict, seed=123):

        parameter_population = []

        self.seed = seed
        np.random.seed(self.seed)

        # Unpack fluxes and concentration into arrays consitent with the
        # compiled functions

        fluxes = [flux_dict[this_reaction.name] for this_reaction in
                  compiled_model.reactions.values()]
        concentrations = np.array([concentration_dict[this_variable] for
                  this_variable in compiled_model.variables.keys()])

        symbolic_concentrations_dict = {Symbol(k):v
                                        for k,v in concentration_dict.items()}

        trials = 0
        while (len(
                parameter_population) < self.parameters.n_samples) or trials > 1e6:
            try:

                parameter_sample = self._sample_saturations_step(compiled_model,
                                                                 symbolic_concentrations_dict,
                                                                 flux_dict)
            except ValueError:
                continue


            # Check stability: real part of all eigenvalues of the jacobian is <= 0

            this_jacobian = compiled_model.jacobian_fun(fluxes, concentrations,
                                                        parameter_sample)
            largest_eigenvalue = eigenvalues(this_jacobian, k=1, which='LR',
                                             return_eigenvectors=False)
            is_stable = largest_eigenvalue <= 0

            compiled_model.logger.info('Model is stable? {} '
                                       '(max real part eigv: {}'.
                                       format(is_stable,largest_eigenvalue))

            if is_stable:
                parameter_population.append(parameter_sample)

            # Count the trials
            trials += 1

        return parameter_population

    def _sample_saturations_step(self, compiled_model, concentration_dict,
                                 flux_dict):
        parameter_sample = {v.symbol: v.value for k,v in compiled_model.parameters.items()}
        # Update the concentrations which are parameters (Boundaries)
        for k,v in concentration_dict.items():
            parameter_sample[k] = v

        # Sample parameters for every reaction
        for this_reaction in compiled_model.reactions.values():


            vmax_param = this_reaction.parameters.vmax_forward

            try:
                keq_param = this_reaction.parameters.k_equilibrium
                this_parameters = {
                    keq_param.symbol: keq_param.value,
                    vmax_param.symbol: 1.0}

            except AttributeError:
                this_parameters = { vmax_param.symbol: 1.0}

            # Get parameters from mechanism
            this_reaction_parameters = this_reaction.parameters

            # Add concentrations
            for k, v in this_reaction_parameters.items():
                try:
                    this_reaction_parameters[k].value = concentration_dict[v.symbol]
                except KeyError:
                    pass

            # Loop over the named tuple
            for this_p_name, this_parameter in this_reaction_parameters.items():
                # Sample a saturation
                # The parameters that have to be sampled are attached to
                # reactants. hence, their .hook attribute shall not be None

                if (this_parameter.hook is not None) \
                   and (this_parameter.value is None):
                    this_saturation = sample()
                    # TODO THIS IS A HOT FIX AND REALLY STUPID REMOVE ASAP
                    # TODO - OK, removed
                    this_reactant = this_parameter.hook.symbol
                    this_concentration = concentration_dict[this_reactant]
                    this_param_symbol = this_parameter.symbol
                    this_parameters[this_param_symbol] = \
                        ((1.0 - this_saturation) * this_concentration) / this_saturation

            # Calculate the effective saturation
            this_net_reaction_rate = this_reaction.mechanism.reaction_rates[
                'v_net']
            this_parameter_subs = concentration_dict.copy()
            this_parameter_subs.update(this_parameters.copy())

            normed_net_reaction_rate = this_net_reaction_rate.evalf(
                subs=this_parameter_subs )

            # If is 0 try more exact evaluation
            if normed_net_reaction_rate == 0:
                normed_net_reaction_rate = float(
                    this_net_reaction_rate.subs(this_parameter_subs))


            if (flux_dict[this_reaction.name] > 0 and
               normed_net_reaction_rate <= 0) \
               or \
               (flux_dict[this_reaction.name] < 0 and
                normed_net_reaction_rate >= 0):
                    # msg = 'Overall saturation for reaction {}' \
                    #       ' is not 0 < {} < 1 '.format(this_reaction.name,
                    #                                    normed_net_reaction_rate)
                    msg = 'Reaction {} operates in opposite direction '.format(this_reaction.name)

                    compiled_model.logger.error(msg)
                    raise ValueError(msg)


            # Calculate the effective VMax
            this_vmax = flux_dict[this_reaction.name] / normed_net_reaction_rate
            this_parameters[vmax_param.symbol] = this_vmax

            # Update the dict with explicit model parameters
            parameter_sample.update(this_parameters)
        return parameter_sample
