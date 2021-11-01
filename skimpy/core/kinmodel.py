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
from skimpy.analysis.ode.utils import make_gamma_fun
from skimpy.analysis.ode.symbolic_jacobian_fun import SymbolicJacobianFunction

from skimpy.analysis.mca.make import make_mca_functions
from skimpy.analysis.mca.prepare import prepare_mca
from skimpy.analysis.mca import *
from skimpy.analysis.mca.volume_ratio_function import VolumeRatioFunction
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
        # Add using add compartments!
        self.compartments = iterable_to_tabdict([])

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

        # Compartment parameters
        for this_comp in self.compartments.values():
            comp_params = TabDict({str(p.symbol): p for p in this_comp.parameters.values()})
            parameters.update(comp_params)

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

        # Compartment parameters
        for this_comp in self.compartments.values():
            comp_params = TabDict({str(p.symbol): p for p in this_comp.parameters.values()})
            parameters.update(comp_params)

        for key,value in value_dict.items():
            parameters[str(key)].value = value

    @property
    def moieties(self):
        ms = []
        if hasattr(self,'conservation_relation'):
            for row in self.conservation_relation.toarray():
                ms.append({self.reactants.iloc(j)[0]: e for j, e in enumerate(row) if not e == 0})
        else:
            raise AttributeError('Model not prepared for MCA')

        return ms

    def add_reaction(self, reaction):
        """
        Adds a SKiMPy reaction to the model

        :param reaction: The reaction to add
        :type reaction: skimpy.core.Reaction
        :return:
        """
        # If the variable name already exists substitute the reactant
        # with the pre-existing variable
        for k,v in reaction.reactants.items():
            if v.name in self.reactants.keys():

                # TODO substitute by itemsetter in reactions.reactants
                # Possible reactant types for the reaction.reactant, all handled by the if/elif/else conditions in order
                # 1. Small molecule  2. Simple activator 3. Simple inhibitor 4. Competitive inhibitor 5. Normal reactant
                # This way, the inhibitors/activators don't get added to reaction.mechanism.reactants and are correctly
                # linked
                if k.startswith('small_molecule'):
                    for this_mod in reaction.modifiers.values():
                        if 'small_molecule' in this_mod.reactants.keys():
                            if this_mod.reactants['small_molecule'].name \
                               is v.name:

                               this_mod.reactants['small_molecule'] = self.reactants[v.name]
                elif k.startswith('activator_'):
                    for this_mod in reaction.modifiers.values():
                        if 'activator' in this_mod.reactants.keys():
                            if this_mod.reactants['activator'].name \
                                    is v.name:
                                this_mod.reactants['activator'] = self.reactants[v.name]
                elif k.startswith('inhibitor_'):
                    for this_mod in reaction.modifiers.values():
                        if 'inhibitor' in this_mod.reactants.keys():
                            if this_mod.reactants['inhibitor'].name \
                                    is v.name:
                                this_mod.reactants['inhibitor'] = self.reactants[v.name]
                elif k.startswith('inhibitor'):
                    for this_inh_name, this_inh in reaction.mechanism.inhibitors.items():
                        if this_inh.name is v.name:
                            reaction.mechanism.inhibitors[this_inh_name] = self.reactants[v.name]
                else:
                    reaction.mechanism.reactants[k] = self.reactants[v.name]


        self.add_to_tabdict(reaction, 'reactions')

    def add_compartment(self, compartment):
        """

        :param compartment:
        :return:
        """
        self.add_to_tabdict(compartment, 'compartments')


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

    def repair(self):

        """
        Link inhibitors and activators to reactants
        FIXME: Any idea to avoid this is dearly welcome
        :return:
        """
        for this_reaction in self.reactions.values():
            this_mechanism = this_reaction.mechanism
            if not this_mechanism.inhibitors is None:
                for this_keys, this_inhibitor in this_mechanism.inhibitors.items():
                    if this_inhibitor.name in self.reactants:
                        this_mechanism.inhibitors[this_keys] = self.reactants[this_inhibitor.name]

            for this_modifier in this_reaction.modifiers.values():
                for this_keys, this_modifier_reactant in this_modifier.reactants.items():
                    if this_modifier_reactant.name in self.reactants:
                        this_modifier.reactants[this_keys] = self.reactants[this_modifier_reactant.name]


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

        if self.compartments:
            self.volume_ratio_func = VolumeRatioFunction(self,
                                                         self.variables,
                                                         self.parameters, )
        else:
            self.volume_ratio_func = None

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


    def compile_jacobian(self, type=NUMERICAL ,sim_type=QSSA, ncpu=1):

        self.sim_type = sim_type

        if not hasattr(self, 'pool'):
            self.pool = Pool(ncpu)

        if type == NUMERICAL:
            self.compile_mca(parameter_list=[], sim_type=sim_type, ncpu=ncpu)

        if type == SYMBOLIC:
            self.compile_ode(sim_type=sim_type, ncpu=ncpu)
            self.jacobian_fun = SymbolicJacobianFunction(self.ode_fun.variables,
                                                         self.ode_fun.expressions,
                                                         self.parameters,
                                                         self.pool)

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

        #Update fixed parameters
        self.ode_fun.get_params()

        # #if parameters are empty try to fetch from model
        # if not self.ode_fun._parameter_values:
        #     self.ode_fun.parameter_values = {v.symbol:v.value
        #                                      for k,v in self.parameters.items()}

        # solve the ode
        solution = self.solver.solve(time_out, ordered_initial_conditions)

        return ODESolution(self, solution)

    def compile_mca(self, parameter_list=[], mca_type=NET, sim_type=QSSA, ncpu=1):
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
                                         sim_type=sim_type,
                                         mca_type=mca_type,
                                        )

                self.independent_elasticity_fun = independent_elasticity_fun
                self.dependent_elasticity_fun = dependent_elasticity_fun
                self.parameter_elasticities_fun = parameter_elasticities_fun

                if mca_type == SPLIT:
                    self.displacement_function = make_gamma_fun(self)
                else:
                    self.displacement_function = None

                # Build functions for stability and control coefficient's
                self.jacobian_fun = JacobianFunction(
                    self.reduced_stoichiometry,
                    self.independent_elasticity_fun,
                    self.dependent_elasticity_fun,
                    self.volume_ratio_func,
                    self.conservation_relation,
                    self.independent_variables_ix,
                    self.dependent_variables_ix)

                self.concentration_control_fun = ConcentrationControlFunction(
                    self,
                    self.reduced_stoichiometry,
                    self.independent_elasticity_fun,
                    self.dependent_elasticity_fun,
                    self.parameter_elasticities_fun,
                    self.volume_ratio_func,
                    self.conservation_relation,
                    self.independent_variables_ix,
                    self.dependent_variables_ix,
                    displacement_function=self.displacement_function,
                    mca_type=mca_type)

                self.flux_control_fun = FluxControlFunction(
                    self,
                    self.reduced_stoichiometry,
                    self.independent_elasticity_fun,
                    self.dependent_elasticity_fun,
                    self.parameter_elasticities_fun,
                    self.conservation_relation,
                    self.independent_variables_ix,
                    self.dependent_variables_ix,
                    self.concentration_control_fun,
                    mca_type=mca_type)
