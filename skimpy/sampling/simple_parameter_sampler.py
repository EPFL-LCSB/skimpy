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

from collections import namedtuple
import numpy as np
from numpy.random import sample
# from scipy.sparse.linalg import eigs as eigenvalues
from scipy.linalg import eigvals as eigenvalues
from sympy import Symbol

import scipy

from skimpy.utils.namespace import *

from skimpy.sampling import ParameterSampler, SaturationParameterFunction, FluxParameterFunction


class SimpleParameterSampler(ParameterSampler):
    """
    A simple parameter sampler that samples stable model parameters
    with respect to a steady state flux and concentration state
    """

    Parameters = namedtuple('Parameters', ['n_samples'])
    # TODO Talk to Pierre / Misko about simple sampler parameters
    # if parameters are not defined put default values
    Parameters.__new__.__defaults__ = (None,) * len(Parameters._fields)

    def sample(self,
               compiled_model,
               flux_dict,
               concentration_dict,
               only_stable=True,
               min_max_eigenvalues=False,
               seed=123,
               bounds_sample=(0, 1),
               max_trials=1e6):

        parameter_population = []
        smallest_eigenvalues = []
        largest_eigenvalues = []

        self.seed = seed
        np.random.seed(self.seed)
        self.bounds_sample = bounds_sample

        # Unpack fluxes and concentration into arrays consitent with the
        # compiled functions

        fluxes = [flux_dict[this_reaction.name] for this_reaction in
                  compiled_model.reactions.values()]
        concentrations = np.array([concentration_dict[this_variable] for
                                   this_variable in
                                   compiled_model.variables.keys()])

        symbolic_concentrations_dict = {Symbol(k): v
                                        for k, v in concentration_dict.items()}

        # Compile functions for calculating Km's and Vmax's
        if not hasattr(compiled_model, 'saturation_parameter_function')\
           or not hasattr(compiled_model, 'flux_parameter_function'):
            self._compile_sampling_functions(
                compiled_model,
                symbolic_concentrations_dict,
                flux_dict)

        # Sample as long as max trials aren't exceeded and we haven't reached
        # target population size
        trials = 0
        while (len(
                parameter_population) < self.parameters.n_samples) and trials < max_trials:

            parameter_sample = self._sample_saturation_step_compiled(
                compiled_model,
                symbolic_concentrations_dict,
                flux_dict)

            # Check stability: real part of all eigenvalues of the jacobian is
            # <= 0
            this_jacobian = compiled_model.jacobian_fun(fluxes, concentrations,
                                                        parameter_sample)

            # largest_eigenvalue = eigenvalues(this_jacobian, k=1, which='LR',
            #                                  return_eigenvectors=False)
            # Test suggests that this is apparently much faster ....

            #this_real_eigenvalues = sorted(np.real(eigenvalues(this_jacobian.todense())))

            # Next level optimization using direclty the lapack function
            lamr, lami, vl, vr, info = scipy.linalg.lapack.dgeev(this_jacobian.todense(),
                                                                 compute_vl=1,
                                                                 compute_vr=0,
                                                                 )

            this_real_eigenvalues = sorted(lamr)

            largest_eigenvalue = this_real_eigenvalues[-1]
            smallest_eigenvalue = this_real_eigenvalues[0]

            is_stable = largest_eigenvalue <= 0

            compiled_model.logger.info('Model is stable? {} '
                                       '(max real part eigv: {})'.
                                       format(is_stable, largest_eigenvalue))

            if is_stable or not only_stable:
                parameter_population.append(parameter_sample)
                largest_eigenvalues.append(largest_eigenvalue)
                smallest_eigenvalues.append(smallest_eigenvalue)

            # Count the trials
            trials += 1

        if min_max_eigenvalues:
            return parameter_population, largest_eigenvalues, smallest_eigenvalues
        else:
            return parameter_population

    # Under construction new sampling with compiled function
    def _compile_sampling_functions(self, model,
                                    concentrations,
                                    fluxes):
        """
        Compiles the function for sampling using cython
        :param model:
        """
        model.saturation_parameter_function = SaturationParameterFunction(model,
                                                                          model.parameters,
                                                                          concentrations)

        model.flux_parameter_function = FluxParameterFunction(model,
                                                              model.parameters,
                                                              concentrations,)

    def _sample_saturation_step_compiled(self,
                                         compiled_model,
                                         concentration_dict,
                                         flux_dict,
                                         parameters_to_resample=None,
                                         fixed_parameters=None):

        """
        Sample one set of saturations using cython complied functions
        :param compiled_model:
        :param concentration_dict:
        :param flux_dict:
        :return:
        """

        parameter_sample = {v.symbol: v.value for k, v in compiled_model.parameters.items()}

        model_parameters = compiled_model.parameters
        # Update the concentrations which are parameters (Boundaries)
        for k, v in concentration_dict.items():
            if str(k) in model_parameters:
                parameter_sample[k] = v

        # Set all vmax/flux parameters to 1.
        # TODO Generalize into Flux and Saturation parameters
        for this_reaction in compiled_model.reactions.values():
            try:
                if this_reaction.enzyme is None:
                    vmax_param = this_reaction.parameters.vmax_forward
                    parameter_sample[vmax_param.symbol] = 1
                else:
                    vmax_param = this_reaction.parameters.kcat_forward
                    parameter_sample[vmax_param.symbol] = 1

            except AttributeError:
                pass

        if not hasattr(compiled_model, 'saturation_parameter_function')\
           or not hasattr(compiled_model, 'flux_parameter_function'):
            raise RuntimeError("Function for sampling not complied")

        # Sample appropriate # of saturations
        if not compiled_model.saturation_parameter_function.sym_saturations:
            n_stats = 0
            saturations = []
        else:
            n_sats = len(compiled_model.saturation_parameter_function.sym_saturations)
            a,b = self.bounds_sample
            saturations = sample(n_sats)*(b-a) + a

        # Calculate the Km's
        compiled_model.saturation_parameter_function(
            saturations,
            parameter_sample,
            concentration_dict,
            parameters_to_resample,
            fixed_parameters
        )

        # Calculate the Vmax's
        compiled_model.flux_parameter_function(
            compiled_model,
            parameter_sample,
            concentration_dict,
            flux_dict
        )

        return parameter_sample
