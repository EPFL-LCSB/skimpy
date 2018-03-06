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
from scipy.integrate import ode
from sympy import sympify

from skimpy.analysis.ode.ode_fun import ODEFunction
from skimpy.analysis.ode.flux_fun import FluxFunction
from skimpy.utils import iterable_to_tabdict
from skimpy.utils.general import join_dicts


def make_ode_fun(kinetic_model, sim_type):

    # Get all variables and expressions (Better solution with types?)
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

    all_expr, all_param = list(zip(*all_data))

    # Flatten all the lists
    flatten_list = lambda this_list: [item for sublist in this_list \
                                      for item in sublist]

    all_rates = flatten_list([these_expressions.keys()
                              for these_expressions in all_expr])
    all_param = join_dicts(all_param)

    # Get unique set of all the variables
    variables = [sympify(x) for x in set(all_rates)]
    variables = iterable_to_tabdict(variables, use_name=False)

    expr = dict.fromkeys(variables, 0.0)

    # Mass balance
    # Sum up all rate expressions
    for this_reaction in all_expr:
        for this_variable_key in this_reaction:
            expr[this_variable_key] += this_reaction[this_variable_key]

    # Apply constraints. Constraints are modifiers that act on
    # expressions
    for this_constraint in kinetic_model.constraints.values():
        this_constraint(expr)

    # Apply boundary conditions. Boundaries are modifiers that act on
    # expressions
    for this_boundary_condition in kinetic_model.boundary_conditions.values():
        this_boundary_condition(expr)

    # Make vector function from expressions
    ode_fun = ODEFunction(variables, expr, all_param)

    return ode_fun, variables


def make_flux_fun(kinetic_model):

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
    # are all the epxressions indexed by the
    flux_fun = FluxFunction(variables, expr, all_param)

    return flux_fun



def get_ode_solver(  ode_fun,
                     solver_type = "vode",
                     reltol = 1e-8,
                     abstol = 1e-8):

    # Initialize the integrator
    ode_solver = ode(ode_fun)
    # Set properties
    ode_solver.set_integrator(  solver_type,
                                method='bdf',
                                atol=abstol,
                                rtol=reltol )

    return ode_solver


def _solve_ode(solver, time_int, initial_concentrations):

    solver.set_initial_value(initial_concentrations, time_int[0])

    t_sol = [time_int[0]]
    y_sol = [initial_concentrations]

    while solver.t <= time_int[1] and solver.successful():
            solver.integrate(time_int[1], step=True)
            t_sol.append(solver.t)
            y_sol.append(solver.y)

    return t_sol,y_sol



