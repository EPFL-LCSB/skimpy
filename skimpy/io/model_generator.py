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
from abc import ABC, abstractmethod
from collections import namedtuple
from skimpy.core import *
from skimpy.utils.namespace import WATER_FORMULA

MetWithStoich = namedtuple('MetWithStoich', ['metabolite', 'stoichiometry'])

class ModelGenerator(ABC):
    def __init__(self,
                 reaction_to_mechanisms=None,
                 reactant_relations=None,
                 small_molecules=None,
                 small_molecule_modifier=None,
                 reactants_to_exclude=None,
                 reaction_groups=None):
        """
        This class defines the rules to build a kinetic models from
        solely stoichiometric data and supplemented data

        :param reaction_to_mechanisms: definition of kinetic mechanisms
        :param reactant_relations: relations of reactants e.g. S1 to P1 and P2
        :param small_molecules: list of small molecules ids

        """
        if reaction_to_mechanisms is not None:
            self.reaction_to_mechanisms = reaction_to_mechanisms
        else:
            self.reaction_to_mechanisms = dict()

        self.reactant_relations = reactant_relations
        self.reaction_groups = reaction_groups
        if reaction_groups is not None:
            self.reactions_to_reaction_groups = {rxn: group for group, list in reaction_groups.items() for rxn in list}
        else:
            self.reactions_to_reaction_groups = None
        """Default values """
        if small_molecules is None:
            self.small_molecules = \
                ['co2', 'nh4', 'o2', 'pi', 'ppi', 'pppi', 'h2o2', 'hco3', 'h2s', 'so3', 'so4']
        else:
            self.small_molecules = small_molecules

        if small_molecule_modifier is None:
            self.small_molecule_modifier = DisplacementSmallMoleculeModifier
        else:
            self.small_molecule_modifier = small_molecule_modifier

        if reactants_to_exclude is None:
            self.reactants_to_exclude = []
        else:
            self.reactants_to_exclude = reactants_to_exclude

        self.dummy_dgo = -10.0


    @abstractmethod
    def import_model(self,model):
        pass

    @abstractmethod
    def import_reaction(self,reaction):
        pass

