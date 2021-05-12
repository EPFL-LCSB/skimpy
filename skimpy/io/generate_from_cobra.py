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
from .model_generator import ModelGenerator, MetWithStoich
from .utils import create_reaction_from_stoich, check_boundary_reaction
from skimpy.utils.general import sanitize_cobra_vars
from skimpy.utils.namespace import WATER_FORMULA


class FromCobra(ModelGenerator):
    """
    Class to generate Kinetic models from cobra
    """
    def __init__(self,
                 **kwargs):

        ModelGenerator.__init__(self,
                                **kwargs
                                )

    def import_model(self, cobra_model):
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
        for this_reaction in cobra_model.reactions:
            if not check_boundary_reaction(this_reaction):
                this_kinetic_reaction = self.import_reaction(cobra_model, this_reaction)
                if this_kinetic_reaction is not None:
                    this_mechanism = this_kinetic_reaction.mechanism
                    parameters[this_kinetic_reaction.name] = this_mechanism.Parameters()
                    skimpy_model.add_reaction(this_kinetic_reaction)

        # Add Boundaries
        for this_reaction in cobra_model.reactions:

            if check_boundary_reaction(this_reaction):
                for this_met in this_reaction.metabolites:
                    # If the metabolite does not correspond to water as water is
                    # omitted from the reactions or if we force the reactant to
                    # be excluded
                    if not (this_met .formula == WATER_FORMULA) \
                       and (this_met .id not in self.reactants_to_exclude):
                        met = sanitize_cobra_vars(this_met.id)
                        this_reactant = skimpy_model.reactants[met]
                        this_const_met = ConstantConcentration(this_reactant)
                        skimpy_model.add_boundary_condition(this_const_met)

        skimpy_model.parametrize_by_reaction(parameters)
        return skimpy_model

    def import_reaction(self, cobra_model, cobra_reaction, name=None, irrev_direction=0):

        if name is None:
            name = cobra_reaction.id

        # Ignore if only water is participating
        is_water = all([met.formula == WATER_FORMULA
                        for met in cobra_reaction.metabolites])

        # Ignore if only reactants participate that will be ignored
        is_ignored = all([met.id in self.reactants_to_exclude
                          for met in cobra_reaction.metabolites])

        if is_ignored or is_water:
            return None

        try:
            reaction_data = self.reaction_to_mechanisms[name]
            skimpy_reaction = create_reaction_from_data(name,
                                                        reaction_data)
        except KeyError:
            met_stoich_dict = {}
            for this_met, this_stoich in cobra_reaction.metabolites.items():
                this_met_id = sanitize_cobra_vars(this_met)
                met_stoich_dict[this_met_id] = MetWithStoich(metabolite=this_met,
                                                             stoichiometry=this_stoich)
            # Get inhibitors from reaction groups
            if self.reactions_to_reaction_groups is not None:
                try:
                    reaction_group = self.reactions_to_reaction_groups[name]
                    reactions_in_group = self.reaction_groups[reaction_group]
                    reactants = set([])
                    products  = set([])
                    this_reactants = set(cobra_reaction.reactants)
                    this_products  = set(cobra_reaction.products)

                    for rxn_id in reactions_in_group:
                        this_cobra_reaction = cobra_model.reactions.get_by_id(rxn_id)
                        reactants.update(this_cobra_reaction.reactants)
                        products.update(this_cobra_reaction.products)

                    inhibitors = products.difference(this_products)\
                        .union(reactants.difference(this_reactants))
                    inhibitors = [i for i in inhibitors]

                except KeyError:
                    inhibitors = None
                    reaction_group = None
            else:
                inhibitors = None
                reaction_group = None

            if irrev_direction > 0:
                # Irreversible in forward direction
                irrev = True
            elif irrev_direction < 0:
                # Irreversible in backward direction
                # invert the assumed forward stoichiometry
                irrev = True
                for this_met_id, metstoich in met_stoich_dict.items():

                    raise NotImplementedError('FIX ME SEE COMMENT')
                    # FIXME: The attribute of the namedtuple cant be changed
                    # find solution for this!!!

                    metstoich.stoichiometry *= -1

            else:
                # Not irreverisble
                irrev = False
                pass

            skimpy_reaction = create_reaction_from_stoich(name,
                                                          met_stoich_dict,
                                                          self,
                                                          inhibitors=inhibitors,
                                                          reaction_group=reaction_group,
                                                          irrev=irrev)

        return skimpy_reaction
