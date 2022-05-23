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

from math import log

from skimpy.core import *
from skimpy.utils.conversions import deltag0_to_keq
from skimpy.utils.general import sanitize_cobra_vars
from skimpy.utils.namespace import WATER_FORMULA
from .model_generator import ModelGenerator, MetWithStoich
from .generate_from_cobra import FromCobra

from .utils import create_reaction_from_stoich, check_boundary_reaction


class FromPyTFA(FromCobra):
    """
    Class to generate Kinetic models from cobra
    """
    def __init__(self,
                 max_revesible_deltag_0 = 100,
                 **kwargs ):

        ModelGenerator.__init__(self,
                                **kwargs
                                )

        self.max_revesible_deltag_0 = max_revesible_deltag_0

    def import_model(self, pytfa_model, pytfa_solution_raw, concentration_scaling_factor = 1.0):
        """
        Function to create a kinetic model from a constraint based model

        :param pytfa_model:
        :param pytfa_solution: a prepresentative solution for the pytfa model solution.raw
        :return: skimpy model
        """

        skimpy_model = KineticModel()

        """Read the reactions """
        # DM_     Boundary reactions
        # By default the metabolites of boundary reactions
        # to be constant concentrations
        parameters = {}
        for this_reaction in pytfa_model.reactions:
            if not check_boundary_reaction(this_reaction):

                k_eq, delta_g0, delta_g = self.get_equlibrium_constant(pytfa_model,
                                                             pytfa_solution_raw,
                                                             this_reaction,
                                                             scaling_factor=concentration_scaling_factor)

                # Check the reversibility of the reaction
                forward = delta_g < 0
                backward = delta_g > 0
                large_delta_g0 = abs(delta_g0) > self.max_revesible_deltag_0
                if forward and large_delta_g0:
                    # Forward
                    irrev_direction = 1
                elif backward and large_delta_g0:
                    # Barckwards
                    irrev_direction = -1
                else:
                    irrev_direction = 0

                this_skimpy_reaction = self.import_reaction(pytfa_model,
                                                            this_reaction,
                                                            name=this_reaction.id,
                                                            irrev_direction=irrev_direction)

                if this_skimpy_reaction is not None:

                    this_mechanism = this_skimpy_reaction.mechanism
                    parameters[this_skimpy_reaction.name] = this_mechanism.Parameters(k_equilibrium=k_eq)

                    skimpy_model.add_reaction(this_skimpy_reaction)

        for this_reaction in pytfa_model.reactions:
            # TODO Check wtf id vs name in pytfa
            if check_boundary_reaction(this_reaction):
                for this_met in this_reaction.metabolites:
                    # If the metabolite does not correspond to water as water is
                    # omitted from the reactions or if we force the reactant to
                    # be excluded
                    if  not (this_met .formula == WATER_FORMULA) \
                        and (this_met .id not in self.reactants_to_exclude):
                        try:
                            this_reactant = skimpy_model.reactants[this_met]
                            this_const_met = ConstantConcentration(this_reactant)
                            skimpy_model.add_boundary_condition(this_const_met)
                        except KeyError:
                            skimpy_model.logger.warning('Metabolite {} in pyTFA model is not part of any reaction and is '
                                                        'omitted in the SKiMpy model '
                                                        .format(this_met,))

        skimpy_model.parametrize_by_reaction(parameters)

        return skimpy_model

    def get_equlibrium_constant(self, pytfa_model, pytfa_solution_data, this_reaction, scaling_factor = 1.0):
        # get delta_Gstd variable name from LC and Delta G
        temp = pytfa_model.TEMPERATURE
        gas_constant = pytfa_model.GAS_CONSTANT
        RT = pytfa_model.RT
        # We here calculate the delta G from
        try:
            scaling_order = sum(this_reaction.metabolites.values() )
            var_delta_g = pytfa_model.delta_g.get_by_id(this_reaction.id).name
            deltag0 = pytfa_solution_data[var_delta_g]
            deltag  = pytfa_solution_data[var_delta_g]
            is_in_model = True
        except KeyError:
            is_in_model = False

        if is_in_model:
            # Calculate the deltaG0 based on the reactants that will
            # be part in the model
            for met, s in pytfa_model.reactions.get_by_id(this_reaction.id).metabolites.items():
                if not (met.formula == WATER_FORMULA) and \
                   met.id not in self.reactants_to_exclude:
                    var_met_lc = pytfa_model.log_concentration.get_by_id(met.id).name
                    met_lc = pytfa_solution_data[var_met_lc]
                    deltag0 -= s * RT * (met_lc + log(scaling_factor))

        else:
            deltag0 = self.dummy_dgo
            deltag  = 0

        k_eq = deltag0_to_keq(deltag0,
                              temp,
                              gas_constant=gas_constant)

        return k_eq, deltag0, deltag







