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
from .model_generator import ModelGenerator
from .utils import create_reaction_from_stoich, check_boundary_reaction
from skimpy.utils.general import sanitize_cobra_vars

class FromCobra(ModelGenerator):
    """
    Class to generate Kinetic models from cobra
    """
    def __init__(self,
                 reaction_to_mechanisms={},
                 reactant_relations={},
                 small_molecules=None,
                 water=None,
                 hydrogen=None,
                 ):
        ModelGenerator.__init__(self,
                                reaction_to_mechanisms=reaction_to_mechanisms,
                                reactant_relations=reactant_relations,
                                small_molecules=small_molecules,
                                water=water,
                                hydrogen=hydrogen)

    def import_model(self,cobra_model):
        """
        Function to create a kinetic model from a constraint based model
        :param cobra_model:
        :return: skimpy model
        """

        skimpy_model = KineticModel()

        """Read the reactions """
        # By default the metabolites of boundary reactions
        # to be constant concentrations
        for this_reaction in cobra_model.reactions:
           if not check_boundary_reaction(this_reaction):
                this_kinetic_reaction = self.import_reaction(this_reaction)
                if this_kinetic_reaction is not None:
                    skimpy_model.add_reaction(this_kinetic_reaction)

        # Add Boundaries
        for this_reaction in cobra_model.reactions:
            if check_boundary_reaction(this_reaction):
                met = sanitize_cobra_vars(this_met.id)

                # If the metabolite does not correspond to water as water is
                # omitted from the reactions
                if not met.startswith("{}_".format(self.water)) \
                and not met.startswith("{}_".format(self.hydrogen)):
                    this_reactant = skimpy_model.reactants[met]
                    this_const_met = ConstantConcentration(this_reactant)
                    skimpy_model.add_boundary_condition(this_const_met)

        return skimpy_model

    def import_reaction(self, cobra_reaction, name=None):

        if name is None:
            name = cobra_reaction.name

        # Ignore if only water is participating
        is_water = all([met.id.startswith("{}_".format(self.water))
                        for met in cobra_reaction.metabolites])

        # Ignore if only hydrogen is participating
        is_hydrogen = all([met.id.startswith("{}_".format(self.hydrogen))
                        for met in cobra_reaction.metabolites])


        if is_hydrogen:
            return None

        try:
            reaction_data = self.reaction_to_mechanisms[name]
            skimpy_reaction = create_reaction_from_data(name,
                                                        reaction_data)
        except KeyError:
            met_stoich_dict = {}
            for k,v in cobra_reaction.metabolites.items():
                k = sanitize_cobra_vars(k.id)
                met_stoich_dict[k] = v

            skimpy_reaction = create_reaction_from_stoich(name,
                                                          met_stoich_dict,
                                                          self)

        return skimpy_reaction




