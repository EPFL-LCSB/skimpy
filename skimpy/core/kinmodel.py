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
from numpy import array, append

from .ode_fun import ODEFunction
from .solution import Solution

from ..utils import TabDict
from ..utils.general import join_dicts

from sympy import sympify

def iterable_to_tabdict(iterable, use_name = True):
    """
    Takes the items from an iterable and puts them in a TabDict, indexed by the
    elements' .name property

    :param iterable:
    :return:
    """
    if iterable is None:
        return TabDict()

    if use_name:
        return TabDict([(x.name, x) for x in iterable])
    else:
        return TabDict([(x.__str__(), x) for x in iterable])

class KineticModel(object):
    # Consult with Pierre about this class!
    # Better use dicts with names! + inherited objects!

    def __init__(self,
                 reactions=None,
                 boundary_conditions=None,
                 constraints=None):
        # initialize the model is stated
        # FIXME Add dictlists from cobra ? or reimplement a similar data structure
        # self.metabolites = metabolites    #List of metabolite objects/ids

        # List of enzyme objects
        self.reactions   = iterable_to_tabdict(reactions)
        self.boundary_conditions  = iterable_to_tabdict(boundary_conditions)
        self.constraints  = iterable_to_tabdict(constraints)
        self.initial_conditions = iterable_to_tabdict([])
        self._modified    = True

    # TODO : Implement
    @property
    def metabolites(self):
        pass

    def add_reaction(self, reaction):
        self.add_to_tabdict(reaction, 'reactions')

    def add_constraint(self, constraint):
        constraint.link(self)
        self.add_to_tabdict(constraint, 'constraints')

    def add_boundary_condition(self, boundary_condition):
        boundary_condition.link(self)
        self.add_to_tabdict(boundary_condition, 'boundary_conditions')

    def add_to_tabdict(self, element, kind):

        the_tabdict = getattr(self, kind)

        # Add an enzyme to the model
        if element.name in the_tabdict:
            error_msg = 'Reaction {} already exists in the model'
            raise(Exception(error_msg.format(element.name)))

        the_tabdict[element.name] = element

        # TODO : Implement with metabolites
        # for this_metabolite in reaction.metabolites:
        #     self.metabolites.append(this_metabolite)
        self._modifed = True


    def parametrize(self,param_dict):
        """
        If has input: apply as dict to reactions in the model by
            reaction.parametrize(args)
        :return:
        """

        for reaction_name, the_params in param_dict.items():
            the_reaction = self.reactions[reaction_name]
            the_reaction.parametrize(the_params)

    @property
    def sim_type(self):
        return self._simtype

    @sim_type.setter
    def sim_type(self, value):
        self._simtype = value
        self._modified = True

    def compile_ode(self,
                    sim_type='QSSA',):

        self.sim_type = sim_type

        # Recompile only if modified
        if self._modified:
            # Compile ode function
            self.make_ode_fun(sim_type)
            self._modified = False
            # Create initial_conditions from variables
            self.initial_conditions = TabDict([(x,0.0) for x in self.variables])

    def solve_ode(self,
                  time_int,
                  solver_type = 'vode',
                  reltol = 1e-8,
                  abstol = 1e-8):

        # Choose a solver
        self.solver = get_ode_solver(self.ode_fun,solver_type,reltol,abstol)

        # Order the initial conditions according to variables
        ordered_initial_conditions = [self.initial_conditions[variable]
                                    for variable in self.variables]

        # solve the ode
        t_sol,y_sol = _solve_ode(self.solver,
                                 time_int,
                                 ordered_initial_conditions)

        return Solution(self,t_sol,y_sol)

    def make_ode_fun(self, sim_type):

        # Gete all variables and expressions (Better solution with types?)
        if sim_type == 'QSSA':
            all_data = [this_reaction.mechanism.get_qssa_rate_expression()  \
                        for this_reaction in self.reactions.values()]

        elif sim_type == 'tQSSA':
            raise(NotImplementedError)
            all_data = [this_reaction.mechanism.get_tqssa_rate_expression() \
                        for this_reaction in self.reactions.values()]

        elif sim_type == 'full':
            all_data = [this_reaction.mechanism.get_full_rate_expression()  \
                        for this_reaction in self.reactions.values()]

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
        for this_constraint in self.constraints.values():
            this_constraint(expr)

        # Apply boundary conditions. Boundaries are modifiers that act on
        # expressions
        for this_boundary_condition in self.boundary_conditions.values():
            this_boundary_condition(expr)


        # Make vector function from expressions
        self.ode_fun = ODEFunction(variables, expr, all_param)
        # Ode variables
        self.variables = variables

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
