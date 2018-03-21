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


def make_mca_functions(kinetic_model,parameter_list,sim_type):
    """ Create the elasticity and flux functions for MCA
    :param kinmodel:
    :param parameter_list:
    :return:
    """

    # Get all variables and expressions (Better solution with types?)
    # TODO This should be a method in KineticModel that stores the expressions
    if sim_type == 'QSSA':
        all_data = [this_reaction.mechanism.get_qssa_rate_expression() \
                    for this_reaction in kinetic_model.reactions.values()]

    elif sim_type == 'tQSSA':
        raise(NotImplementedError)
        all_data = [this_reaction.mechanism.get_tqssa_rate_expression() \
                    for this_reaction in kinetic_model.reactions.values()]

    elif sim_type == 'full':
        all_data = [this_reaction.mechanism.get_full_rate_expression() \
                    for this_reaction in kinetic_model.reactions.values()]


    # Get flux expressions for the net
    all_flux_expressions = [this_reaction.mechanism.reaction_rates['v_net'] \
                           for this_reaction in kinetic_model.reactions]

    all_expr, all_parameters = list(zip(*all_data))

    # Flatten all the lists
    flatten_list = lambda this_list: [item for sublist in this_list \
                                      for item in sublist]

    all_rates = flatten_list([these_expressions.keys()
                              for these_expressions in all_expr])

    # Sort into an ordered list
    all_parameters = [sympify(x) for x in join_dicts(all_parameters).keys()]
    all_parameters = iterable_to_tabdict(all_parameters, use_name=False)

    # Get unique set of all the variables
    all_variables = [sympify(x) for x in set(all_rates)]
    all_variables = iterable_to_tabdict(all_variables, use_name=False)

    #TODO handle depedent variables (i.e. concentrations)
    relative_weights = []
    absolute_weights = []

    #######

    all_independent_variables = []

    all_dependent_variables = []

    #parameter elasticity function
    parameter_elasticities_fun = make_elasticity_fun(all_flux_expressions,
                                                     parameter_list,
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
           parameter_elasticities_fun,\
           relative_weights,\
           absolute_weights, \
           all_variables, \
           all_parameters



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

