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
from skimpy.core import *

class ModelGenerator(ABC):
    def __init__(self,
                 reaction_to_mechanisms=None,
                 reactant_relations=None,
                 small_molecules=None,
                 small_molecule_modifier=None,
                 water=None,
                 hydrogen=None):
        """
        This class defines the rules to build a kinetic models from
        solely stoichiometric data and supplemented data
        :param reaction_to_mechanisms: definition of kinetic mechanisms
        :param reactant_relations: relations of reactants e.g. S1 to P1 and P2
        :param small_molecules: list of small molecules ids
        :param water: id of the water molecule
        """
        self.reaction_to_mechanisms = reaction_to_mechanisms
        self.reactant_relations = reactant_relations


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

        if water is None:
            self.water = 'h2o'
        else:
            self.water = water

        if hydrogen is None:
            self.hydrogen = 'h'
        else:
            self.hydrogen = hydrogen


    @abstractmethod
    def import_model(self,model):
        pass

    @abstractmethod
    def import_reaction(self,reaction):
        pass

