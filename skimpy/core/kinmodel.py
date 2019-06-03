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

from scikits.odes import ode
from skimpy.analysis.ode.utils import make_ode_fun
from skimpy.analysis.mca.make import make_mca_functions
from skimpy.analysis.mca.prepare import prepare_mca
from skimpy.analysis.mca import *
from ..utils.logger import get_bistream_logger
from .solution import ODESolution

from ..utils import TabDict, iterable_to_tabdict
from ..utils.namespace import *

from multiprocessing import Pool

class KineticModel(object):
    """
    This class contains the kinetic model as described by reaction and
    boundary conditions and constratins.

    :param :
    :return:
    """

    def __init__(self,
                 reactions=None,
                 boundary_conditions=None,
                 constraints=None,
                 name='Unnamed'):
        self.name = name
        self.reactions = iterable_to_tabdict(reactions)
        self.boundary_conditions = iterable_to_tabdict(boundary_conditions)
        self.constraints = iterable_to_tabdict(constraints)
        self.initial_conditions = iterable_to_tabdict([])
        self.logger = get_bistream_logger(name)
        self._simtype = None
        self._modified = True
        self._recompiled = False

        #self.parameters = TabDict()

    @property
    def reactants(self):
        reactants = TabDict([])
        for this_reaction in self.reactions.values():
            this_rectants = TabDict([(v.name,v) for v in this_reaction.reactants.values()])
            reactants.update(this_rectants)
        return reactants

    @property
    def parameters(self):
        parameters = TabDict([])
        for this_reaction in self.reactions.values():
            reaction_params = TabDict({str(p.symbol): p for p in this_reaction.parameters.values()})
            parameters.update(reaction_params)
        return parameters

    @parameters.setter
    def parameters(self,value_dict):
        """

        :param value_dict:
        :return: Nothing
        """
        parameters = TabDict([])
        for this_reaction in self.reactions.values():
            reaction_params = TabDict({str(p.symbol): p for p in this_reaction.parameters.values()})
            parameters.update(reaction_params)

        for key,value in value_dict.items():
            parameters[str(key)].value = value


    def add_reaction(self, reaction):
        """
        Adds a SKiMPy reaction to the model

        :param reaction: The reaction to add
        :type reaction: skimpy.core.Reaction
        :return:
        """
        # If the variable name already exists substitute
        # with the variable
        for k,v in reaction.reactants.items():
            if v.name in self.reactants.keys():

                # TODO substitute by itemsetter in reactions.reactants
                # UGLYYYYYYY
                if k.startswith('small_molecule'):
                    for this_mod in reaction.modifiers.values():
                        if this_mod.reactants['small_molecule'].name \
                           is v.name:

                           this_mod.reactants['small_molecule'] = self.reactants[v.name]
                else:
                    reaction.mechanism.reactants[k] = self.reactants[v.name]


        self.add_to_tabdict(reaction, 'reactions')

    def add_constraint(self, constraint):
        constraint.link(self)
        self.add_to_tabdict(constraint, 'constraints')

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

    def parametrize_by_reaction(self, param_dict):
        """
        If has input: apply as dict to reactions in the model by
            reaction.parametrize(args)
        :return:
        """
        for reaction_name, the_params in param_dict.items():
            the_reaction = self.reactions[reaction_name]
            the_reaction.parametrize(the_params)

        #self.update()

    def parametrize(self, param_dict):
        raise NotImplemented('We have to do some thinking OK')
        pass


    @property
    def sim_type(self):
        return self._simtype

    @sim_type.setter
    def sim_type(self, value):
        self._simtype = value
        self._modified = True

    def prepare(self, mca=True, ode=True, **kwargs):
        """
        Model preparation for different analysis types. The preparation is done before the compiling step
        to be able to curate the model in between

        :param mca:
        :param ode:
        :return:
        """
        self.variables = TabDict([(k, v.symbol) for k, v in self.reactants.items()])

        if mca:
            reduced_stoichometriy,\
            conservation_relation, \
            independent_variables_ix, \
            dependent_variables_ix = prepare_mca(kinetic_model=self, **kwargs)

            self.conservation_relation = conservation_relation
            self.reduced_stoichiometry = reduced_stoichometriy
            self.dependent_variables_ix = dependent_variables_ix
            self.independent_variables_ix = independent_variables_ix

            self.dependent_reactants = TabDict([self.reactants.iloc(ix)
                                        for ix in self.dependent_variables_ix])

        if ode:
            pass

    def compile_ode(self, sim_type=QSSA, ncpu=1,):

        # For security
        # self.update()

        self.sim_type = sim_type

        if not hasattr(self, 'pool'):
            self.pool = Pool(ncpu)

        # Recompile only if modified or simulation
        if self._modified or self.sim_type != sim_type:
            # Compile ode function
            ode_fun, variables = make_ode_fun(self, sim_type, pool=self.pool)
            # TODO define the init properly
            self.ode_fun = ode_fun
            self.variables = variables

            self._modified = False
            self._recompiled = True
            # Create initial_conditions from variables
            old_initial_conditions = self.initial_conditions
            self.initial_conditions = TabDict([(x,0.0) for x in self.variables])
            # If data was stored previously in the initial conditions, recover it (needed for
            # serialization)
            self.initial_conditions.update(old_initial_conditions)

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
        if not hasattr(self, 'solver')\
           or self._recompiled:
            self.solver = ode(solver_type, self.ode_fun, **kwargs)
            self._recompiled = False

        # Order the initial conditions according to variables
        ordered_initial_conditions = [self.initial_conditions[variable]
                                      for variable in self.variables]

        # #if parameters are empty try to fetch from model
        # if not self.ode_fun._parameter_values:
        #     self.ode_fun.parameter_values = {v.symbol:v.value
        #                                      for k,v in self.parameters.items()}

        # solve the ode
        solution = self.solver.solve(time_out, ordered_initial_conditions)

        return ODESolution(self, solution)

    def compile_mca(self, parameter_list=[], sim_type=QSSA, ncpu=1):
            """
            Compile MCA expressions: elasticities, jacobian
            and control coeffcients
            """
            if not hasattr(self, 'pool'):
                self.pool = Pool(ncpu)

            # Recompile only if modified or simulation
            if self._modified or self.sim_type != sim_type:

                # Get the model expressions
                independent_elasticity_fun, \
                dependent_elasticity_fun, \
                parameter_elasticities_fun, \
                    = make_mca_functions(self,
                                         parameter_list,
                                         sim_type=sim_type
                                        )

                self.independent_elasticity_fun = independent_elasticity_fun
                self.dependent_elasticity_fun = dependent_elasticity_fun
                self.parameter_elasticities_fun = parameter_elasticities_fun


                # Build functions for stability and control coefficient's
                self.jacobian_fun = JacobianFunction(
                    self.reduced_stoichiometry,
                    self.independent_elasticity_fun,
                    self.dependent_elasticity_fun,
                    self.conservation_relation,
                    self.independent_variables_ix,
                    self.dependent_variables_ix)

                self.concentration_control_fun = ConcentrationControlFunction(
                    self,
                    self.reduced_stoichiometry,
                    self.independent_elasticity_fun,
                    self.dependent_elasticity_fun,
                    self.parameter_elasticities_fun,
                    self.conservation_relation,
                    self.independent_variables_ix,
                    self.dependent_variables_ix)

                self.flux_control_fun = FluxControlFunction(
                    self,
                    self.reduced_stoichiometry,
                    self.independent_elasticity_fun,
                    self.dependent_elasticity_fun,
                    self.parameter_elasticities_fun,
                    self.conservation_relation,
                    self.independent_variables_ix,
                    self.dependent_variables_ix,
                    self.concentration_control_fun)
