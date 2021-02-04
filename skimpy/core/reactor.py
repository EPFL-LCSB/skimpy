# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2020 Laboratory of Computational Systems Biotechnology (LCSB),
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

from copy import copy
from abc import ABC
from scikits.odes import ode

from skimpy.core.itemsets import Reactant
from skimpy.utils.namespace import *
from skimpy.utils.medium import get_medium

from skimpy.utils import TabDict, iterable_to_tabdict
from skimpy.analysis.ode.utils import get_expressions_from_model, make_expressions
from skimpy.analysis.ode.ode_fun import ODEFunction

from skimpy.core.solution import ODESolution

from multiprocessing.pool import Pool

from sympy import exp, Symbol

class Reactor(ABC):
    """

    """
    def __init__(self,
                 models,
                 biomass_reactions,
                 biomass_scaling,
                 boundary_conditions=None,
                 extracellular_compartment='e',
                 custom_variables=[],
                 ):
        """

        :param models:
        :param biomass_reactions: Dict, TabDict
        :param boundary_conditions:
        :param extracellular_compartment:
        """

        # Medium of the reactor
        # Copy models, Generate unique model annotation and remove all boundary conditions
        models = [copy(m) for m in models]
        for the_model in models:
            the_model.boundary_conditions.clear()

        medium = []
        for model in models:
            medium += get_medium(model, extracellular_compartment=extracellular_compartment)
        self.extracellular_compartment = extracellular_compartment

        self.medium = iterable_to_tabdict(medium)

        self.biomass_reactions = biomass_reactions
        self.biomass_scaling = biomass_scaling

        # biomass variables
        biomass_variables = [ Reactant('biomass_{}'.format(model.name), model=model)
                             for model in models]
        self.biomass_variables = iterable_to_tabdict(biomass_variables)

        #Generate unique model annotation and remove all boundary conditions
        for the_model in models:
            reactants = the_model.reactants
            for the_reactant in reactants.values():
                if the_reactant.name not in self.medium:
                    the_reactant.name = the_model.name + '_' + the_reactant.name
                    the_reactant._generate_symbol()

            parameters = the_model.parameters
            for the_parameter in parameters.values():
                the_parameter.name = the_model.name + '_' + the_parameter.name
                the_parameter._generate_symbol()

            for the_reaction in the_model.reactions.values():
                the_reaction.name = the_model.name + '_' + the_reaction.name

        self.models = iterable_to_tabdict(models)

        self.boundary_conditions = iterable_to_tabdict(boundary_conditions)
        self.initial_conditions = iterable_to_tabdict([])

        self._modified = True
        self.custom_variables = iterable_to_tabdict(custom_variables)

    @property
    def variables(self):
        variables = copy(self.biomass_variables)
        variables.update(self.medium)
        for this_model in self.models.values():
            this_rectants = TabDict([(v.name, v) for v in this_model.reactants.values()
                                     if v.name not in self.medium])
            variables.update(this_rectants)
            # Add the biomass
        variables.update(self.custom_variables)
        return variables

    @property
    def parameters(self):
        parameters = TabDict([])
        for this_model in self.models.values():
            this_parameters = TabDict([(str(v.symbol), v) for v in this_model.parameters.values()])
            parameters.update(this_parameters )
        return parameters


    def parametrize(self, value_dict, model_name):
        """
        Set the parameters of a specified model
        :param value_dict:
        :param model_name:
        :return:
        """
        model = self.models[model_name]

        parameters = TabDict([])
        for this_reaction in model.reactions.values():
            reaction_params = TabDict({str(p.symbol): p for p in this_reaction.parameters.values()})
            parameters.update(reaction_params)

            # Compartment parameters
        for this_comp in model.compartments.values():
            comp_params = TabDict({str(p.symbol): p for p in this_comp.parameters.values()})
            parameters.update(comp_params)

        for key, value in value_dict.items():
            if model_name+'_'+str(key) in parameters.keys():
                parameters[model_name+'_'+str(key)].value = value


    def add_boundary_condition(self, boundary_condition):
        """
        Enforces a boundary condition (e.g. a constant concentration) on the
        kinetic model

        :param boundary_condition: the boundary condition to enforce
        :type boundary_condition: skimpy.core.BoundaryCondition
        :return:
        """
        boundary_condition.link(self)
        self.add_to_tabdict(boundary_condition, 'boundary_conditions')

        # HOTFIX to be done better
        boundary_condition.reactant.type = VARIABLE

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
        self._modified = True

    def compile_ode(self, sim_type=QSSA, ncpu=1, add_dilution=False, custom_ode_update=None):
        """

        :param sim_type:
        :param ncpu:
        :return:
        """
        self.sim_type = sim_type
        if not hasattr(self, 'pool'):
            self.pool = Pool(ncpu)

        # Recompile only if modified or simulation
        if self._modified or self.sim_type != sim_type:
            # Compile ode function
            ode_fun, variables = make_reactor_ode_fun(self, sim_type, pool=self.pool,
                                                      add_dilution=add_dilution,
                                                      custom_ode_update=custom_ode_update)
            # TODO define the init properly
            self.ode_fun = ode_fun
            self._modified = False
            self._recompiled = True
            # Create initial_conditions from variables
            old_initial_conditions = self.initial_conditions
            self.initial_conditions = TabDict([(x,0.0) for x in self.variables])
            # If data was stored previously in the initial conditions, recover it (needed for
            # serialization)
            self.initial_conditions.update(old_initial_conditions)

    def initialize(self, value_dict, model_name):
        """
        Set the initial conditions for a single strain (excluding medium)
        :param value_dict:
        :param model_name:
        :return:
        """
        for key,value in value_dict.items():
            if not key in self.medium:
                if model_name+'_'+str(key) in self.initial_conditions:
                    self.initial_conditions[model_name+'_'+str(key)] = value


    def solve_ode(self, time_out, solver_type='cvode', **kwargs):
        """
        The solver types are from ::scikits.odes::, and can be found at
        <https://scikits-odes.readthedocs.io/en/latest/solvers.html>`_.
        :param time_out: The times at which the solution is evaluated
        :type time_out:  list(float) or similar
        :param solver_type: must be among ['cvode','ida','dopri5','dop853']
        :param kwargs:
        :return:
        """
        extra_options = {'old_api': False}
        kwargs.update(extra_options)

        # Choose a solver
        if not hasattr(self, 'solver') \
                or self._recompiled:
            self.solver = ode(solver_type, self.ode_fun, **kwargs)
            self._recompiled = False

        # Order the initial conditions according to variables
        ordered_initial_conditions = [self.initial_conditions[variable]
                                      for variable in self.variables]
        #Update fixed parameters
        self.ode_fun.get_params()
        solution = self.solver.solve(time_out, ordered_initial_conditions)

        return ODESolution(self, solution)


def make_reactor_ode_fun(reactor, sim_type, pool=None, add_dilution=False,
                         custom_ode_update=None):
    """
    This function generates the symbolic expressions for a reactor model

    :param reactor:
    :param sim_type:
    :param pool:
    :return:
    """

    expression_list = []
    parameters_list = []

    # Flatten all the lists
    flatten_list = lambda this_list: [item for sublist in this_list \
                                      for item in sublist]
    medium_symbols = [m.symbol for m in reactor.medium.values()]

    model_biomass_var = {b.model:b for b in reactor.biomass_variables.values()}

    for model in reactor.models.values():
        medium_com = reactor.extracellular_compartment
        volume_scaling_medium = model.compartments[medium_com].parameters.cell_volume.symbol / \
                                model.compartments[medium_com].parameters.volume.symbol

        biomass_symbol = model_biomass_var[model].symbol*volume_scaling_medium

        all_data = get_expressions_from_model(model, sim_type,
                                              medium_symbols=medium_symbols,
                                              biomass_symbol=biomass_symbol)
        # get the dxdt expressions 
        all_expr_model, _, all_parameters_model = list(zip(*all_data))
        parameters_list.extend(all_parameters_model)
        expression_list.extend(all_expr_model)

    parameters_list = flatten_list(parameters_list)
    parameters_list = list(set(parameters_list))
    parameters_list = iterable_to_tabdict(parameters_list, use_name=False)

    variables = TabDict([(k,v.symbol) for k,v in reactor.variables.items()])

    volume_ratios = TabDict([])

    for model in reactor.models.values():
        # Volumes
        if model.compartments:
            this_volume_ratios = TabDict([(k,v.compartment.parameters.cell_volume.symbol/
                                           v.compartment.parameters.volume.symbol)
                                          if k not in reactor.medium else (k, 1.0)
                                          for k,v in model.reactants.items()
                                          ])
            for comp in model.compartments.values():
                this_comp_parameters = {str(v.symbol):v.symbol for v in comp.parameters.values() }
                parameters_list.update( this_comp_parameters )
        else:
            this_volume_ratios = {}
        volume_ratios.update(this_volume_ratios)

    for k,biomass in reactor.biomass_variables.items():
        volume_ratios[k] = 1.0

    for k,var in reactor.custom_variables.items():
            volume_ratios[k] = 1.0

    # Make default ODE expressions
    expr = make_expressions(variables,expression_list, volume_ratios=volume_ratios, pool=pool)

    if add_dilution:
        # Dilution for intracellular metabolites
        for model in reactor.models.values():
            growth_reaction = reactor.biomass_reactions[model.name]
            growth_rate_expression = growth_reaction.mechanism.reaction_rates['v_net']
            biomass_scaling = reactor.biomass_scaling[model.name]

            for r in model.reactants.values():
                if not r.name in reactor.medium:
                    expr[r.symbol] -= r.symbol*growth_rate_expression/biomass_scaling

    # Add growth expressions
    for biomass in reactor.biomass_variables.values():
        growth_reaction = reactor.biomass_reactions[biomass.model.name]
        growth_rate_expression = growth_reaction.mechanism.reaction_rates['v_net']
        biomass_scaling = reactor.biomass_scaling[biomass.model.name]
        expr[biomass.symbol] += biomass.symbol*growth_rate_expression/biomass_scaling

    # Apply constraints. Constraints are modifiers that act on
    # expressions
    for model in reactor.models.values():
        for this_constraint in model.constraints.values():
            this_constraint(expr)

    # Apply boundary conditions. Boundaries are modifiers that act on
    # expressions

    for this_boundary_condition in reactor.boundary_conditions.values():
        this_boundary_condition(expr)

    # Make vector function from expressions
    ode_fun = ODEFunction(reactor, variables, expr, parameters_list, pool=pool,
                          custom_ode_update=custom_ode_update)

    return ode_fun, variables

