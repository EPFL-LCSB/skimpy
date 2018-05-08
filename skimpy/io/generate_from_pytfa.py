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

from skimpy.core import *
from skimpy.utils.conversions import deltag0_to_keq
from .model_generator import ModelGenerator
from .generate_from_cobra import FromCobra

from .utils import create_reaction_from_stoich


class FromPyTFA(FromCobra):
    """
    Class to generate Kinetic models from cobra
    """
    def __init__(self,
                 reaction_to_mechanisms={},
                 reactant_relations={},
                 small_molecules=None,
                 water=None):
        ModelGenerator.__init__(self,
                                reaction_to_mechanisms=reaction_to_mechanisms,
                                reactant_relations=reactant_relations,
                                small_molecules=small_molecules,
                                water=water)

    def import_model(self,pytfa_mode):
        """
        Function to create a kinetic model from a constraint based model
        :param cobra_model:
        :return: skimpy model
        """

        skimpy_model = KineticModel()

        """Read the reactions """
        # DM_     Boundary reactions
        # By default the metabolites of boundary reactions
        # to be constant concentrations
        parameters = {}
        for this_reaction in pytfa_mode.reactions:

            if this_reaction.id.startswith("DM_"):
                for this_met in this_reaction.metabolites:
                    this_const_met = ConstantConcentration(this_met.name)
                    skimpy_model.add_boundary_condition(this_const_met)
            else:

                this_skimpy_reaction = self.import_reaction(this_reaction,name = this_reaction.id)

                #get delta G
                deltag0 = getattr(pytfa_mode.delta_g,
                                  this_reaction.id).variable.primal

                temp = pytfa_mode.TEMPERATURE
                gas_constant = pytfa_mode.GAS_CONSTANT
                k_eq = deltag0_to_keq(deltag0,
                                      temp,
                                      gas_constant=gas_constant)

                this_mechanism=this_skimpy_reaction.mechanism
                parameters[this_skimpy_reaction.name]=this_mechanism.Parameters(k_equilibrium=k_eq)
                skimpy_model.add_reaction(this_skimpy_reaction)

        skimpy_model.parametrize(parameters)

        return skimpy_model







