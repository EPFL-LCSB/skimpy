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
from sympy import simplify

from skimpy.analysis.ode.ode_fun import ODEFunction
from skimpy.analysis.ode.flux_fun import FluxFunction
from skimpy.analysis.ode.gamma_fun import GammaFunction

from skimpy.utils import iterable_to_tabdict, TabDict
from skimpy.utils.namespace import *
from skimpy.utils.general import join_dicts


def make_ode_fun(kinetic_model, sim_type, pool=None, custom_ode_update=None):
    """

    :param kinetic_model:
    :param sim_type:
    :return:
    """
    all_data = get_expressions_from_model(kinetic_model, sim_type)

    # get expressions for dxdt
    all_expr, _, all_parameters = list(zip(*all_data))

    # Flatten all the lists
    flatten_list = lambda this_list: [item for sublist in this_list \
                                      for item in sublist]

    all_parameters = flatten_list(all_parameters)
    all_parameters = list(set(all_parameters))
    all_parameters = iterable_to_tabdict(all_parameters, use_name=False)


    # Better since this is implemented now
    reactant_items = kinetic_model.reactants.items()
    variables = TabDict([(k,v.symbol) for k,v in reactant_items])

    #Compartments # CHECK IF THIS ONLY IS TRUE IF ITS NOT EMPTY
    if kinetic_model.compartments:
        #TODO Throw error if no cell reference compartment is given

        volume_ratios = TabDict([(k,v.compartment.parameters.cell_volume.symbol/
                            v.compartment.parameters.volume.symbol )
                           for k,v in kinetic_model.reactants.items()])
        for comp in kinetic_model.compartments.values():
            this_comp_parameters = {str(v.symbol):v.symbol for v in comp.parameters.values() }
            all_parameters.update( this_comp_parameters )
    else:
        volume_ratios = None

    expr = make_expressions(variables,all_expr, volume_ratios=volume_ratios ,pool=pool)

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
    ode_fun = ODEFunction(kinetic_model, variables, expr, all_parameters, pool=pool,
                          custom_ode_update=custom_ode_update)

    return ode_fun, variables


def make_flux_fun(kinetic_model, sim_type):
    """

    :param kinetic_model:
    :param sim_type:
    :return:
    """
    all_data = get_expressions_from_model(kinetic_model, sim_type)

    # Get flux expressions
    _, all_expr, all_parameters = list(zip(*all_data))

    reactions = kinetic_model.reactions.keys()

    # Flatten all the lists
    flatten_list = lambda this_list: [item for sublist in this_list \
                                      for item in sublist]

    if sim_type == ELEMENTARY:
        expr = [[(r+'_'+er, ex) for er, ex in e.items()] for r, e in zip(reactions, all_expr)]
        expr = TabDict(flatten_list(expr))
    else:
        expr = TabDict([(r, e) for r, e in zip(reactions, all_expr)])


    all_parameters = flatten_list(all_parameters)
    all_parameters = list(set(all_parameters))
    all_parameters = iterable_to_tabdict(all_parameters, use_name=False)

    # Better since this is implemented now
    reactant_items = kinetic_model.reactants.items()
    variables = TabDict([(k,v.symbol) for k,v in reactant_items])

    # Make vector function from expressions in this case all_expressions
    # are all the expressions indexed by the
    flux_fun = FluxFunction(variables, expr, all_parameters)
    flux_fun._parameter_values = {v:p.value for v,p in kinetic_model.parameters.items()}

    return flux_fun


def make_gamma_fun(kinetic_model):
    """
    Return a function that calculates the thermodynamic displacement for
    all the reactions in a model
    :param kinetic_model:
    :return:
    """
    reactions = kinetic_model.reactions.values()
    all_parameters = []
    expr = TabDict([])
    for r in reactions:
        try:
            keq = r.parameters.k_equilibrium.symbol
            expr[r.name] = 1/keq

            # Here we need to get all variables and parameters
            reactant_stoichiometry = TabDict([])
            for k, v in r.mechanism.reactant_stoichiometry.items():
                this_reactant = r.mechanism.reactants[k]
                if this_reactant in reactant_stoichiometry.keys():
                    reactant_stoichiometry[this_reactant] += v
                else:
                    reactant_stoichiometry[this_reactant] = v

            for this_modifier in r.modifiers.values():
                for k,v in this_modifier.reactants.items():
                    if this_modifier.reactant_stoichiometry[k] != 0:
                        reactant_stoichiometry[v] = this_modifier.reactant_stoichiometry[k]

            for v, s in reactant_stoichiometry.items():
                expr[r.name] *= v.symbol**s
                if v.type == PARAMETER:
                    all_parameters.append(v.symbol)

            all_parameters.append(keq)

        except AttributeError:
            pass

    all_parameters = iterable_to_tabdict(all_parameters, use_name=False)

    # Better since this is implemented now
    reactant_items = kinetic_model.reactants.items()
    variables = TabDict([(k,v.symbol) for k,v in reactant_items])

    # Make vector function from expressions in this case all_expressions
    # are all the expressions indexed by the
    gamma_fun = GammaFunction(variables, expr, all_parameters)
    gamma_fun._parameter_values = {v:p.value for v,p in kinetic_model.parameters.items()}

    return gamma_fun



def make_expressions(variables, all_flux_expr, volume_ratios=None,pool=None):

    if pool is None:
        expr = dict.fromkeys(variables.values(), 0.0)

        for this_reaction in all_flux_expr:
            for this_variable_key in this_reaction:
                try:
                    expr[this_variable_key] += this_reaction[this_variable_key]
                except KeyError:
                    pass

    else:

        inputs = [(v,all_flux_expr) for v in variables.values()]

        list_expressions = pool.map(make_expresson_single_var, inputs)

        expr = join_dicts(list_expressions)

    #Add compartment volumes
    if not volume_ratios is None:
        for k,v in variables.items():
            volume_ratio = volume_ratios[k]
            # Mutiply massbalance for each metabolite by volume ratio
            expr[v] = volume_ratio*expr[v]


    return expr

def make_expresson_single_var(input):
    var, all_flux_expr = input

    this_expr = dict()
    this_expr[var] = 0.0

    for this_reaction in all_flux_expr:
        for this_variable_key in this_reaction:
            if this_variable_key == var:
                this_expr[this_variable_key] += this_reaction[this_variable_key]

    return this_expr


def get_expressions_from_model(kinetic_model, sim_type,
                               medium_symbols=None,
                               biomass_symbol=None):
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

            # For reactor building
            if not medium_symbols is None:
                vars_in_medium = [v for v in dxdt if v in medium_symbols]
                for v in vars_in_medium:
                    dxdt[v] = dxdt[v]*biomass_symbol

            all_data.append((dxdt, flux, parameters))

    elif sim_type == TQSSA:
        raise(NotImplementedError)

    elif sim_type == ELEMENTARY:
        all_data = []
        #TODO Modifiers sould be applicable for all simulation types
        for this_reaction in kinetic_model.reactions.values():
            this_reaction.mechanism.get_full_rate_expression()

            all_data.append((this_reaction.mechanism.expressions,
                             this_reaction.mechanism.reaction_rates,
                             this_reaction.mechanism.expression_parameters)
                            )
    else:
        raise(ValueError('Simulation type not recognized: {}'.format(sim_type)))

    return all_data