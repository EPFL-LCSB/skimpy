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

from sympy import diff, sympify, simplify

from .elasticity_fun import ElasticityFunction
from ..ode.flux_fun import FluxFunction


def make_mca_functions(kinmodel,mca_type = 'vmax'):
    """"Create the elasticity and flux functions for MCA"""

    # Get flux expressions for the net
    all_flux_expressions = [this_reaction.mechanism.reaction_rates['v_net'] \
                           for this_reaction in kinmodel.reactions]

    all_variables = []

    all_parameters = []

    #TODO handle depedent variables (i.e. concentrations)
    relative_weights = []
    absolute_weights = []

    #######

    all_independent_variables = []

    all_dependent_variables = []

    if type == 'vmax':
        # Get all vmax parameters
        this_parameters = []

    elif type == 'kms':
        raise (NotImplementedError)
    else:
        raise(NotImplementedError)

    #parameter elasticitiesfunction
    parameter_elasticities_fun = make_elasticity_fun(all_flux_expressions,
                                                     this_parameters,
                                                     all_variables,
                                                     all_parameters)

    #concentration elasticity functions
    indepdendent_elasticity_fun = make_elasticity_fun(all_flux_expressions,
                                                      all_independent_variables,
                                                      all_variables,
                                                      all_parameters)

    depdendent_elasticity_fun = make_elasticity_fun(all_flux_expressions,
                                                    all_dependent_variables,
                                                    all_variables,
                                                    all_parameters)


    return indepdendent_elasticity_fun, \
           depdendent_elasticity_fun, \
           parameter_elasticities_fun,
           relative_weights,
           absolute_weights


def make_elasticity_fun(expressions,respective_variables ,variables, parameters):
    """
    Create an ElasticityFunction with elasticity = dlog(expression)/dlog(respective_variable)
    :param expressions  tab_dict of expressions (e.g. forward and backward fluxes)
    :param variables    list of variables as string (e.g. concentrations or parameters)
    
    """

    # Get the derivative of expression x vs variable y
    elasticity_expressions = {}
    column = 0
    row = 0
    for this_expression in expressions.values():
        for this_variable in respective_variables:
            row += 1
            column += 1
            this_elasticity = get_dlogx_dlogy(this_expression, this_variable)

            if this_elasticity != 0:
                elasticity_expressions[(row, column)] = this_elasticity

    # Create the elasticity function
    elasticity_fun = ElasticityFunction(elasticity_expressions, variables, parameters)

    return elasticity_fun


def get_dlogx_dlogy(sympy_expression, string_variable):
    """
    Calc d_log_x/d_log_y = y/x*dx/dy
    """
    variable = sympify(string_variable)
    partial_derivative = diff(sympy_expression, variable)

    expression = simplify(partial_derivative / sympy_expression * variable)

    return

