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

from skimpy.analysis.mca.elasticity_fun import ElasticityFunction
from skimpy.analysis.mca.utils import get_dlogx_dlogy
from skimpy.utils import iterable_to_tabdict
from skimpy.utils.general import join_dicts
from skimpy.utils.namespace import QSSA, PARAMETER, TQSSA, ELEMENTARY, NET, SPLIT

from skimpy.utils import TabDict, iterable_to_tabdict


def make_mca_functions(kinetic_model,parameter_list,sim_type, mca_type=NET):
    """ Create the elasticity and flux functions for MCA
    :param kinmodel:
    :param parameter_list:
    :return:
    """

    # Get all variables and expressions (Better solution with types?)
    # TODO This should be a method in KineticModel that stores the expressions
    if sim_type == QSSA:
        all_data = []
        # TODO Modifiers should be applicable for all simulation types
        for this_reaction in kinetic_model.reactions.values():
            this_reaction.mechanism.get_qssa_rate_expression()
            # Update rate expressions
            for this_mod in this_reaction.modifiers.values():
                this_mod(this_reaction.mechanism.reaction_rates)
            this_reaction.mechanism.update_qssa_rate_expression()

            # Add modifier expressions
            for this_mod in this_reaction.modifiers.values():
                # Get parameters from modifiers
                for p_type, parameter in this_mod.parameters.items():
                    mod_sym = parameter.symbol
                    this_reaction.mechanism.expression_parameters.update([mod_sym])

                for r_type, reactant in this_mod.reactants.items():
                    # Add massbalances for modfier reactants if as non-zero stoich
                    if this_mod.reactant_stoichiometry[r_type] == 0:
                        continue

                    mod_sym = reactant.symbol
                    flux = this_reaction.mechanism.reaction_rates['v_net']
                    flux_expression = flux * this_mod.reactant_stoichiometry[r_type]
                    this_reaction.mechanism.expressions[mod_sym] = flux_expression

                    # Add small molecule parameters if they are
                    if reactant.type == PARAMETER:
                        this_reaction.mechanism.expression_parameters.update([mod_sym])

            flux = this_reaction.mechanism.reaction_rates['v_net']
            dxdt = this_reaction.mechanism.expressions
            parameters = this_reaction.mechanism.expression_parameters

            all_data.append((this_reaction.mechanism.expressions,
                         this_reaction.mechanism.expression_parameters))

    elif sim_type == TQSSA:
        raise(NotImplementedError)

    elif sim_type == ELEMENTARY:
        raise(NotImplementedError)

    else:
        raise ArgumentError('{} is not recognized as a simulation type'.format(sim_type))


    # Get flux expressions for the net
    if mca_type == NET:
        all_flux_expressions = [this_reaction.mechanism.reaction_rates['v_net'] \
                               for this_reaction in kinetic_model.reactions.values()]
    elif mca_type == SPLIT:
        all_flux_expressions_fwd = [this_reaction.mechanism.reaction_rates['v_fwd'] \
                                for this_reaction in kinetic_model.reactions.values()]
        all_flux_expressions_bwd = [this_reaction.mechanism.reaction_rates['v_bwd'] \
                                    for this_reaction in kinetic_model.reactions.values()]
        all_flux_expressions = all_flux_expressions_fwd + all_flux_expressions_bwd

    all_expr, all_parameters = list(zip(*all_data))

    # Flatten all the lists
    flatten_list = lambda this_list: [item for sublist in this_list \
                                      for item in sublist]

    all_rates = flatten_list([these_expressions.keys()
                              for these_expressions in all_expr])

    # Sort into an ordered list
    all_parameters = flatten_list(all_parameters)
    all_parameters = list(set(all_parameters))
    all_parameters = iterable_to_tabdict(all_parameters, use_name=False)

    all_variables = TabDict([(k, v.symbol) for k, v in kinetic_model.reactants.items()])

    all_independent_variables = TabDict([(k,v) for e,(k,v) in enumerate(all_variables.items())
                                         if e in kinetic_model.independent_variables_ix])
    all_dependent_variables = TabDict([(k,v) for e,(k,v) in enumerate(all_variables.items())
                                       if e in kinetic_model.dependent_variables_ix])


    #parameter elasticity function
    if parameter_list:
        parameter_elasticities_fun = make_elasticity_fun(all_flux_expressions,
                                                         parameter_list,
                                                         all_variables,
                                                         all_parameters,
                                                         kinetic_model.pool
                                                         )
    else:
        parameter_elasticities_fun = None

    #concentration elasticity functions
    independent_elasticity_fun = make_elasticity_fun(all_flux_expressions,
                                                      all_independent_variables,
                                                      all_variables,
                                                      all_parameters,
                                                      kinetic_model.pool
                                                     )

    if all_dependent_variables:
        dependent_elasticity_fun = make_elasticity_fun(all_flux_expressions,
                                                        all_dependent_variables,
                                                        all_variables,
                                                        all_parameters,
                                                        kinetic_model.pool
                                                       )
    else:
        dependent_elasticity_fun = None

    return independent_elasticity_fun, dependent_elasticity_fun, parameter_elasticities_fun


def make_elasticity_fun(expressions, respective_variables, variables, parameters, pool=None):
    """
    Create an ElasticityFunction with elasticity = dlog(expression)/dlog(respective_variable)
    :param expressions  tab_dict of expressions (e.g. forward and backward fluxes)
    :param variables    list of variables as string (e.g. concentrations or parameters)

    """
    if pool is None:
        elasticity_fun = make_elasticity_fun_single_cpu(expressions,
                                                       respective_variables,
                                                       variables,
                                                       parameters)
    else:
        elasticity_fun = make_elasticity_fun_multicore(expressions,
                                                      respective_variables,
                                                      variables,
                                                      parameters,
                                                      pool)


    return elasticity_fun


def make_elasticity_fun_single_cpu(expressions,respective_variables ,variables, parameters):
    # Get the derivative of expression x vs variable y
    elasticity_expressions = {}

    for row, this_expression in enumerate(expressions):

        for column,this_variable in enumerate(respective_variables.values()):

            if this_variable in this_expression.free_symbols:
                this_elasticity = get_dlogx_dlogy(this_expression, this_variable)
                elasticity_expressions[(row, column)] = this_elasticity

    # Shape of the matrix
    shape = (len(expressions), len(respective_variables))

    # Create the elasticity function
    elasticity_fun = ElasticityFunction(elasticity_expressions,
                                        respective_variables,
                                        variables,
                                        parameters,
                                        shape, )
    return elasticity_fun


def make_elasticity_fun_multicore(expressions,respective_variables ,variables, parameters, pool):
    # Get the derivative of expression x vs variable y

    inputs = [(i,e,respective_variables) for i,e in enumerate(expressions)]

    all_row_slices = pool.map(make_elasticity_single_row,inputs)

    elasticity_expressions = join_dicts(all_row_slices)

    # Shape of the matrix
    shape = (len(expressions), len(respective_variables))

    # Create the elasticity function
    elasticity_fun = ElasticityFunction(elasticity_expressions,
                                        respective_variables,
                                        variables,
                                        parameters,
                                        shape,
                                        pool=pool)
    return elasticity_fun


def make_elasticity_single_row(input):
    """
    Halter function to compute a full row of the elasticity matrix
    :param input: input tuple of row, expression and all respective vars
    :return:
    """
    elasticity_expressions_row_slice = {}
    this_row, this_expression, respective_variables = input

    for column, this_variable in enumerate(respective_variables.values()):

        if this_variable in this_expression.free_symbols:
            this_elasticity = get_dlogx_dlogy(this_expression, this_variable)
            elasticity_expressions_row_slice[(this_row, column)] = this_elasticity

    return elasticity_expressions_row_slice