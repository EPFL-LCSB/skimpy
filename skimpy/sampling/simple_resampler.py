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
from scipy.linalg import eigvals as eigenvalues
from sympy import Symbol

from skimpy.utils.namespace import *

from skimpy.sampling import SimpleParameterSampler


class SimpleResampler(SimpleParameterSampler):
    """
    A parameter sampler that tries to resample parameters that are not included
    in the given fixed_saturation_population.

    The maximum number of trials to get a stable model is implemented
    differently than `SimpleParameterSampler`. In `SimpleParameterSampler` the
    maximum number of trials is the maximum # of sampling attemps to get a stable
    model when `.sample()` is called. Here, maximum number of trials is defined
    per parameter vector in `fixed_saturation_population`. If a certain parameter
    vector exceeds the maximum # of trials to get a stable model, the method
    moves to the next parameter vector in `fixed_saturation_population`

    Used for performing Global Sensitivity Analysis
    """

    def sample(self,
               compiled_model,
               flux_dict,
               concentration_dict,
               fixed_saturation_population,
               min_max_eigenvalues=False,
               seed=123):

        parameter_population = []
        smallest_eigenvalues = []
        largest_eigenvalues = []

        self.seed = seed
        np.random.seed(self.seed)

        # Unpack fluxes and concentration into arrays consitent with the
        # compiled functions

        fluxes = [flux_dict[this_reaction.name] for this_reaction in
                  compiled_model.reactions.values()]
        concentrations = np.array([concentration_dict[this_variable] for
                                   this_variable in
                                   compiled_model.variables.keys()])

        symbolic_concentrations_dict = {Symbol(k): v
                                        for k, v in concentration_dict.items()}

        trials = 0

        # Compile functions
        self._compile_sampling_functions(
            compiled_model,
            symbolic_concentrations_dict,
            flux_dict)

        # Try to re-sample for each supplied parameter vector
        for this_saturation in fixed_saturation_population:

            # only try 1e4 times to get a stable model, else skip to next
            trials = 0
            while trials < 1e4:
                # try get a stable model
                parameter_sample = self._sample_saturation_step_compiled(
                    compiled_model,
                    symbolic_concentrations_dict,
                    flux_dict)

                # combine the resampled parameters with those from the
                # fixed_saturation_population:
                # 1 - unpack this_saturation to kinetic parameters
                # 2 - merge the this_parameter with parameter_sample

                # Check stability: real part of all eigenvalues of the jacobian
                # is <= 0
                this_jacobian = compiled_model.jacobian_fun(fluxes,
                                                            concentrations,
                                                            parameter_sample)

                # largest_eigenvalue = eigenvalues(this_jacobian, k=1,
                # which='LR', return_eigenvectors=False)
                # Test suggests that this is apparently much faster ....
                this_real_eigenvalues = np.real(sorted(
                    eigenvalues(this_jacobian.todense())))

                largest_eigenvalue = this_real_eigenvalues[-1]
                smallest_eigenvalue = this_real_eigenvalues[0]

                is_stable = largest_eigenvalue <= 0

                compiled_model.logger.info('Model is stable? {} '
                                           '(max real part eigv: {})'.
                                           format(is_stable, largest_eigenvalue))
                if is_stable:
                    parameter_population.append(parameter_sample)
                    largest_eigenvalues.append(largest_eigenvalue)
                    smallest_eigenvalues.append(smallest_eigenvalue)
                    continue
                trials += 1

        if min_max_eigenvalues:
            return parameter_population, largest_eigenvalues, smallest_eigenvalues
        else:
            return parameter_population
