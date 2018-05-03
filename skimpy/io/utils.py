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

from skimpy.mechanisms import *


def create_reaction_from_stoich(name,
                                met_stoich_dict,
                                model_generator):

    water = model_generator.water
    small_molecules = model_generator.small_molecules
    reactant_relations = model_generator.reactant_relations

    this_reaction_small_molecules = {}
    this_reaction_reactatns = {}

    """ 
    Split into small molecules and reactants
    and omit water 
    """

    for this_met, stoich in met_stoich_dict.items():
        this_name = this_met.name
        is_small_molecule = any([this_name.startswith(s) for s in small_molecules])

        if not this_met.name.startswith(water + '_') \
           and is_small_molecule:
            this_reaction_small_molecules[this_met] = stoich

        elif not this_met.name.startswith(water + '_') \
             and not is_small_molecule:
            this_reaction_reactatns[this_met] = stoich

    # Determine mechanism
    mechanism = guess_mechanism(this_reaction_reactatns)
    # Determine reactant order
    reactants = guess_reactant_order(mechanism,this_reaction_reactatns,reactant_relations)

    skimpy_reaction = Reaction(name=name,
                                reactants=reactants,
                                mechanism=mechanism)

    # Add small molecules modifiers
    small_molecule_modifier = model_generator.small_molecule_modifier
    for this_sm in small_molecules:
        small_molecule_mod = small_molecule_modifier(kinetic_reaction, this_sm.name)
        skimpy_reaction.modifiers.append(small_molecule_mod)

    return skimpy_reaction


def create_reaction_from_data(name,
                              reaction_data):
    """

    :param reaction_data:
    :return:
    """
    # TODO determine a data format and implement a parser to save reaction data
    # TODO E.G: GLYC: RandBiBiRevMM, S1 = A S2 = B, P1 = C , P2 = D ......
    raise NotImplmentendError


def guess_mechanism(reactants):
    # 1) Reversible Michaelis Menten kinetics
    # 2) Convenience kinetics
    mechanism = check_rev_michaelis_menten(reactants)
    if mechanism is not None:
        return mechanism

    else:
        return None


def check_rev_michaelis_menten(reactants):
    stoich = set([i for i in this_reaction_reactants.values()])
    if stoich is set([-1,-1,1,1]):
        return RandBiBiReversibleMichaelisMenten

    elif stoich is set([-1,1]):
        return ReversibleMichaelisMenten

    else:
        return None

