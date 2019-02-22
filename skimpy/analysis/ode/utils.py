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
from sympy import sympify

from skimpy.analysis.ode.ode_fun import ODEFunction
from skimpy.analysis.ode.flux_fun import FluxFunction
from skimpy.utils import iterable_to_tabdict, TabDict
from skimpy.utils.namespace import *


def make_ode_fun(kinetic_model, sim_type, ncpu=1):
    """

    :param kinetic_model:
    :param sim_type:
    :return:
    """
    sim_type = sim_type.lower()
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
                small_mol = this_mod.reactants['small_molecule']
                sm = small_mol.symbol
                flux = this_reaction.mechanism.reaction_rates['v_net']
                flux_expression_sm = flux * this_mod.reactant_stoichiometry
                this_reaction.mechanism.expressions[sm] = flux_expression_sm
                # Add small molecule parameters if they are
                if small_mol.type == PARAMETER:
                    this_reaction.mechanism.expression_parameters.update([small_mol.symbol])

            all_data.append((this_reaction.mechanism.expressions,
                             this_reaction.mechanism.expression_parameters))


    elif sim_type == TQSSA:
        raise(NotImplementedError)


    elif sim_type == ELEMENTARY:
        all_data = []
        #TODO Modifiers sould be applicable for all simulation types
        for this_reaction in kinetic_model.reactions.values():
            this_reaction.mechanism.get_full_rate_expression()

            all_data.append((this_reaction.mechanism.expressions,
                             this_reaction.mechanism.expression_parameters)
                            )
    else:
        raise(ValueError('Simulation type not recognized: {}'.format(sim_type)))

    all_expr, all_parameters = list(zip(*all_data))

    # Flatten all the lists
    flatten_list = lambda this_list: [item for sublist in this_list \
                                      for item in sublist]

    all_rates = flatten_list([these_expressions.keys()
                              for these_expressions in all_expr])

    all_parameters = flatten_list(all_parameters)
    all_parameters = list(set(all_parameters))
    all_parameters = iterable_to_tabdict(all_parameters, use_name=False)

    # Get unique set of all the variables
    # variables = [sympify(x) for x in set(all_rates)]
    # variables = iterable_to_tabdict(variables, use_name=False)

    # Better since this is implemented now
    variables = TabDict([(k,v.symbol) for k,v in kinetic_model.reactants.items()])

    #TODO If the model has mojeties only account for the independent vars

    expr = dict.fromkeys(variables.values(), 0.0)

    # Mass balance
    # Sum up all rate expressions
    # Todo Add here the parralel option interm of two functions
    for this_reaction in all_expr:
        for this_variable_key in this_reaction:
            try:
                expr[this_variable_key] += this_reaction[this_variable_key]
            except KeyError:
                pass

    # Apply constraints. Constraints are modifiers that act on
    # expressions
    for this_constraint in kinetic_model.constraints.values():
        this_constraint(expr)

    # NEW: Boundary conditions are now handled as parameters
    # Apply boundary conditions. Boundaries are modifiers that act on
    # expressions
    # for this_boundary_condition in kinetic_model.boundary_conditions.values():
    #     this_boundary_condition(expr)

    # Make vector function from expressions
    ode_fun = ODEFunction(kinetic_model, variables, expr, all_parameters)

    return ode_fun, variables


def make_flux_fun(kinetic_model):
    """

    :param kinetic_model:
    :type kinetic_model: skimpy.core.KineticModel
    :return:
    """

    # Get all variables and expressions (Better solution with types?)
    all_rate_expr = [this_reaction.mechanism.reaction_rates \
                    for this_reaction in kinetic_model.reactions.values()]

    all_param = kinetic_model.ode_fun.parameters


    # Flatten all the lists
    flatten_list = lambda this_list: [item for sublist in this_list \
                                      for item in sublist]

    variables = kinetic_model.ode_fun.variables
    expr = {}

    # Mass balance
    # Sum up all rate expressions
    for this_reaction in all_rate_expr:
        for this_rate_key in this_reaction:
            expr[this_rate_key] = this_reaction[this_rate_key]

    # Make vector function from expressions in this case all_expressions
    # are all the expressions indexed by the
    flux_fun = FluxFunction(variables, expr, all_param)
    flux_fun._parameter_values = kinetic_model.ode_fun.parameter_values

    return flux_fun




