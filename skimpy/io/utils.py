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
from operator import itemgetter
import re

def create_reaction_from_stoich(name,
                                met_stoich_dict,
                                model_generator):

    water = model_generator.water
    hydrogen = model_generator.hydrogen
    small_molecules = model_generator.small_molecules
    reactant_relations = model_generator.reactant_relations

    this_reaction_small_molecules = {}
    this_reaction_reactants = {}

    """ 
    Split into small molecules and reactants
    and omit water 
    """

    for this_met, stoich in met_stoich_dict.items():
        # TODO the detection needs to be better !!!
        is_small_molecule = any([this_met.startswith(s) for s in small_molecules])
        # TODO the detection needs to be better !!!
        if not this_met.startswith("{}_".format(water)) \
           and not this_met.startswith("{}_".format(hydrogen)) \
           and is_small_molecule:

            this_reaction_small_molecules[this_met] = stoich

        elif not this_met.startswith("{}_".format(water)) \
             and not this_met.startswith("{}_".format(hydrogen)) \
             and not is_small_molecule:

            this_reaction_reactants[this_met] = stoich

    # TODO this currently catches transport of small molecules
    # create proper transports
    if not this_reaction_reactants:
        this_reaction_reactants = this_reaction_small_molecules
        this_reaction_small_molecules = {}

    # If there are no valid reactants
    # do not consider the reactions
    if not this_reaction_reactants:
        return None

    # Determine mechanism
    mechanism = guess_mechanism(this_reaction_reactants)
    # Determine reactant order
    reactants = make_reactant_set(mechanism,
                                  this_reaction_reactants,
                                  reactant_relations)

    skimpy_reaction = Reaction(name=name,
                                reactants=reactants,
                                mechanism=mechanism)

    # Add small molecules modifiers
    small_molecule_modifier = model_generator.small_molecule_modifier
    for this_sm, stoich in this_reaction_small_molecules.items():
        small_molecule_mod = small_molecule_modifier(this_sm, stoich)
        skimpy_reaction.modifiers[small_molecule_mod.name] = small_molecule_mod

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
    mechanism = check_rev_michaelis_menten(reactants)
    if mechanism is not None:
        return mechanism
    # Else convineince or irrevv MM
    else:
        stoich = [i for i in reactants.values()]
        stoich.sort(reverse=True)
        if max(stoich) > 3 or len(stoich) > 5:
            return make_irrev_m_n_michaelis_menten(stoich)
        else:
            return make_convenience(stoich)

def check_boundary_reaction(cobra_reaction):
    stoich = [i for i in cobra_reaction.metabolites.values()]
    if all([i < 0for i in stoich] ) or all([i > 0for i in stoich]):
        return True
    else:
        return False


def check_rev_michaelis_menten(reactants):
    stoich = [i for i in reactants.values()]
    stoich.sort()

    if stoich == [-1,-1,1,1]:
        return RandBiBiReversibleMichaelisMenten

    if stoich == [-1,1]:
        return ReversibleMichaelisMenten

    return None


def make_reactant_set(mechanism,
                      reactants,
                      reactant_relations):

    reactants_sorted = TabDict(sorted(reactants.items(),
                                      key=itemgetter(1),
                                      reverse=True))

    reactant_list_temp = sorted(mechanism.reactant_stoichiometry.items(),
                           key=lambda x: (x[1],x[0]), reverse=True)

    reactant_list = [v for v,k in reactant_list_temp]

    # ToDo Sort the reactants to match the reactant list
    # using the reactant_relations

    reactants_resorted = {}
    (mets, stoichiometries) = zip(*reactants_sorted.items())
    reactants_sorted = zip(mets,stoichiometries,reactant_list)

    for met, stoich, react in reactants_sorted:
        reactants_resorted[react] = met

    return mechanism.Reactants(**reactants_resorted)

