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
from numpy.random import sample

from sympy import symbols,Symbol
from skimpy.utils.compile_sympy import make_cython_function

class SaturationParameterFunction():
    """
    A class used in the process of sampling to calculate Km's. Provided with a
    model, creates `self.__call__` function using `Cython` to calculate Km's
    given (sampled) sigmas

    :param model:
    :param parameters: the parameters of the model. Parameters with a `.hook`
    field and an empty `.value` will be sampled
    :param concentrations:
    """
    def __init__(self,model,parameters,concentrations):


        # Gather concentrations and saturation_parameters from the inputs.
        # saturation_parameters are those that will be sampled
        self.sym_concentrations = [c for c in concentrations]
        self.saturation_parameters = [v for k,v in parameters.items()
                                      if (v.hook is not None)
                                      and (v.value is None) ] # TODO: should this be removed?

        if not self.saturation_parameters:
            # If there are no saturation parameters in the model dont compile
            # the function (for example, elementary reaction model)
            self.expressions = None
            self.sym_saturations = None
            self.function = None
        else:
            # Create a cython function that calculates Km's given sigmas and
            # concentrations. First, collect the needed algebraic expressions:
            sym_saturations = []
            expressions = []
            for p in self.saturation_parameters:
                this_sat_symbol = Symbol("sigma_"+str(p.symbol))
                sym_saturations.append(this_sat_symbol )
                expressions.append((1-this_sat_symbol)*p.hook.symbol/this_sat_symbol)

            self.expressions = expressions
            self.sym_saturations = sym_saturations

            # Create the cython function
            sym_vars = sym_saturations + self.sym_concentrations
            self.function = make_cython_function(sym_vars, expressions, simplify=False, pool=model.pool)


    def __call__(self, saturations, parameters, concentrations, parameters_to_resample,
                 fixed_parameters):

        # Transform the sample to bounds accroding to the bounds of the
        # parameters respective to their concentrations
        lower_saturations = []
        upper_saturations = []

        if self.function is None:
            pass
        else:
            for p in self.saturation_parameters:
                # The lower bound of the parameter fixes the upper bound on the
                # concentration and vice versa
                the_lower_bound_saturation = 0.0 if p._upper_bound is None \
                    else concentrations[p.hook.symbol] / \
                         (p._upper_bound + concentrations[p.hook.symbol])

                the_upper_bound_saturation = 1.0 if p._lower_bound is None \
                    else concentrations[p.hook.symbol] / \
                         (p._lower_bound + concentrations[p.hook.symbol])

                lower_saturations.append(the_lower_bound_saturation)
                upper_saturations.append(the_upper_bound_saturation)

            _lower_saturations = np.array(lower_saturations)
            _upper_saturations = np.array(upper_saturations)

            # Scale according to lower/upper bounds. `saturations` are in [0,1]
            _saturations = _lower_saturations + saturations * (upper_saturations - _lower_saturations)

            # Get the numerical values of the concentrations
            _concentrations = np.array([concentrations[c] for c in self.sym_concentrations])

            # Calculate Km's and put result in saturation_parameter_values
            input = np.concatenate((_saturations,_concentrations))
            saturation_parameter_values = np.zeros(len(self.saturation_parameters))
            self.function(input,saturation_parameter_values)

            # Assigning saturation parameters
            if not parameters_to_resample:
                for p,v in zip(self.saturation_parameters, saturation_parameter_values):
                    parameters[p.symbol] = v
            else:
                # Only assign sampled parameters that are in `parameters_to_resample`. Use
                # the value of `fixed_parameters` for other parameters
                for c, p in enumerate(self.saturation_parameters):
                    if p in parameters_to_resample:
                        parameters[p.symbol] = saturation_parameter_values[c]
                    else:
                        parameters[p.symbol] = fixed_parameters[p.symbol]
