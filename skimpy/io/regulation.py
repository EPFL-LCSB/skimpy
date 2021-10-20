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
from skimpy.mechanisms import *
from skimpy.mechanisms.generalized_reversible_hill_n_n_h1_with_inhibition import *
from skimpy.mechanisms.convenience_with_inihibition import *
from skimpy.core.compartments import Compartment
from skimpy.utils.general import get_all_reactants


def load_enzyme_regulation(kmodel, df_regulations_all):
    """
    Test function for loading regulations into an existing kmodel

    :param kmodel: Existing kinetic model
    :param df_regulations_all: a pandas df that has the following columns: reaction_id | regulator | regulation
                                reaction_id must correspond with the names in the kinetic model
                                reactor is the name of a metabolite, with the compartment after the _
                                regulation is either 'inhibition', 'activation', or 'competitive inhibition'
    :return: new_kmodel : A new kinetic model with regulations
    """


    new_kmodel = KineticModel(name=kmodel.name)
    #
    unmodified_rxns = [kmodel.reactions[r] for r in kmodel.reactions if
                       r not in list(df_regulations_all['reaction_id'])]
    modified_rxns = [kmodel.reactions[r] for r in kmodel.reactions if r in list(df_regulations_all['reaction_id'])]

    # Rebuild all unmodified reactions as is
    for this_rxn in unmodified_rxns:
        TheMechanism = this_rxn.mechanism.__class__
        string_reactants = {k:v.name for k,v in this_rxn.mechanism.reactants.items()}
        the_reactants = TheMechanism.Reactants(**string_reactants)

        the_enzyme = this_rxn.enzyme

        new_rxn = Reaction(name=this_rxn.name,
                           mechanism=TheMechanism,
                           reactants=the_reactants,
                           enzyme=the_enzyme,
                           )

        # Add kinetic modifiers - note this has been changed to once again not just copy the old modifier but create
        # it from scratch - Currently, this only covers existing small_molecule modifiers
        # TODO: extend this to all possible modifiers including activation and inhibition
        modifiers = this_rxn.modifiers
        for the_modifier in modifiers.values():
            TheModifier = the_modifier.__class__
            new_modifier = TheModifier(small_molecule=the_modifier.reactants['small_molecule'].name,
                                       mechanism_stoichiometry = the_modifier.reactant_stoichiometry,
                                       reaction=new_rxn)
            new_rxn.modifiers[new_modifier.name] = new_modifier

        new_kmodel.add_reaction(new_rxn)

    # Rebuild modified reactions on a case by case basis
    for this_rxn in modified_rxns:

        rxn_mechanism = str(this_rxn.mechanism)  # TODO: This is a big string - make it nicer
        df_regulations = df_regulations_all[df_regulations_all['reaction_id'] == this_rxn.name]

        # First sort out all competitive inhibitions
        df_competitive = df_regulations[df_regulations['regulation'] == 'competitive inhibition']
        if len(df_competitive) > 0:
            list_competitive_inhibitors = list(df_competitive['regulator'])

            # Get the existing reactant stoichiometry
            reactant_stoichiometry = [v for k, v in this_rxn.mechanism.reactant_stoichiometry.items()]
            inhibitor_stoichiometry = np.ones(len(list_competitive_inhibitors))

            # If the mechanism is Gen Hill or Michaelis Menten then we use Gen Rev Hill w inhibition as mechanism
            if ('ReversibleMichaelisMenten' in rxn_mechanism) or ('H1Generalized' in rxn_mechanism):
                mechanism = make_generalized_reversible_hill_n_n_h1_with_inhibition(reactant_stoichiometry,
                                                                                    inhibitor_stoichiometry)

            else:  # Otherwise, use convenience kinetics if it is uni-bi hill or convenience TODO: expand this list!
                mechanism = make_convenience_with_inhibition(reactant_stoichiometry, inhibitor_stoichiometry)

            # For Rev MM, only 1 substrate / product pair --> needs handling before converting to Gen Hill
            if 'ReversibleMichaelisMenten' in rxn_mechanism:
                reactant_key_value = {k + '1': v.name for k, v in this_rxn.mechanism.reactants.items()}
            elif 'UniBi' in rxn_mechanism:
                reactant_key_value = TabDict([(k + '1', v.name) if k.startswith('substrate') else (k,v.name )
                                              for k, v in this_rxn.mechanism.reactants.items()])
            elif 'BiUni' in rxn_mechanism:
                reactant_key_value = TabDict([(k + '1', v.name) if k.startswith('product') else (k,v.name )
                                              for k, v in this_rxn.mechanism.reactants.items()])
            else:
                reactant_key_value = {k: v.name for k, v in this_rxn.mechanism.reactants.items()}

            inhibitor_key_value = {'inhibitor{}'.format(i + 1): e
                                   for i, e in enumerate(list_competitive_inhibitors)}

            reactants = mechanism.Reactants(**reactant_key_value)
            inhibitors = mechanism.Inhibitors(**inhibitor_key_value)

            parameters = mechanism.Parameters()

            new_rxn = Reaction(name=this_rxn.name,
                               mechanism=mechanism,
                               reactants=reactants,
                               inhibitors=inhibitors,
                               )

            # Add existing kinetic modifiers - note this has been changed to once again not just copy the old modifier but create
            # it from scratch - Currently, this only covers existing small_molecule modifiers
            # TODO: extend this to all possible modifiers including activation and inhibition
            modifiers = this_rxn.modifiers
            for the_modifier in modifiers.values():
                TheModifier = the_modifier.__class__
                new_modifier = TheModifier(small_molecule=the_modifier.reactants['small_molecule'].name,
                                           mechanism_stoichiometry=the_modifier.reactant_stoichiometry,
                                           reaction=new_rxn)
                new_rxn.modifiers[new_modifier.name] = new_modifier


        # If there are no competitive regulations, just recreate the reaction first since the other modifications are applied on top
        else:
            TheMechanism = this_rxn.mechanism.__class__
            string_reactants = {k: v.name for k, v in this_rxn.mechanism.reactants.items()}
            the_reactants = TheMechanism.Reactants(**string_reactants)
            the_enzyme = this_rxn.enzyme

            new_rxn = Reaction(name=this_rxn.name,
                               mechanism=TheMechanism,
                               reactants=the_reactants,
                               enzyme=the_enzyme,
                               )
            # Add kinetic modifiers
            modifiers = this_rxn.modifiers
            for the_modifier in modifiers:
                new_rxn.modifiers[the_modifier] = modifiers[the_modifier]

        # Next go through the rest of the inhibitions/activation
        df_other_regulation = df_regulations[df_regulations['regulation'] != 'competitive inhibition']

        for ix_ in df_other_regulation.index:  # TODO: please write this in a better way..
            this_regulation = df_other_regulation.loc[ix_]['regulation']
            this_regulator = df_other_regulation.loc[ix_]['regulator']

            if this_regulation == 'inhibition':
                new_rxn.modifiers[this_regulator] = InhibitionModifier(this_regulator, reaction=new_rxn)
            elif this_regulation == 'activation':
                new_rxn.modifiers[this_regulator] = ActivationModifier(this_regulator, reaction=new_rxn)

        # Finally add the new reaction to rxns
        new_kmodel.add_reaction(new_rxn)

    # Note DW: This is currently still needed to avoid having 2 objects for the same reactant
    new_kmodel.repair()

    # Resbuild compartments
    for the_comp in kmodel.compartments:
        new_comp = Compartment(the_comp)
        new_kmodel.add_compartment(new_comp)


    # Populate the kinmodel.parameters TabDict
    parameter_init_dict = dict()
    for rxn_obj in new_kmodel.reactions.values():
        # initalize empty param list
        parameter_init_dict[rxn_obj.name] = rxn_obj.mechanism.Parameters()

    new_kmodel.parametrize_by_reaction(parameter_init_dict)

    # Boundary Conditions
    for the_bc in kmodel.boundary_conditions.values():
        TheBoundaryCondition = the_bc.__class__

        # We have the following specifically for 'h_e' which is already a parameter at this stage instead reactant
        # Other ec metabolites are currently reactants. Once they get converted to BCs, they become parameters
        try:
            reactant = new_kmodel.reactants[the_bc.reactant.name]
        except KeyError:
            reactant = new_kmodel.parameters[the_bc.reactant.name]

        # NOTE: Check how to properly get other things for other BCs than CC
        new_bc = TheBoundaryCondition(reactant)
        new_kmodel.add_boundary_condition(new_bc)

        # Do not forget to add the value of the BC! Note BN: it has to be .value unlike for the load_yaml
        reactant.value = kmodel.parameters[str(reactant.symbol)].value

    # Update the parameters with all pre-existing parameter values
    this_params = new_kmodel.parameters
    for parameter, symbol in kmodel.parameters.items():
        try:
            this_params[parameter].value = symbol.value
        except KeyError:
            pass

    # Initial conditions
    # we can copy this since it is of type TabDict((str,float) , )
    new_kmodel.initial_conditions = kmodel.initial_conditions.copy()

    # Adding compartments
    # Note DW: This is needed because not all reactants are variables
    # e.g. extracellular mets with BC == const > thus they are parameters
    # Note BN: This has been moved here from above since now all the parameters have been populated as well
    for the_met in get_all_reactants(kmodel).values():
        try:
            new_met = new_kmodel.reactants[the_met.name]
        except KeyError:
            new_met = new_kmodel.parameters[the_met.name]

        if not the_met.compartment is None:
            comp = new_kmodel.compartments[the_met.compartment.name]
            new_met.compartment = comp

    # If the model the model has computed moieties we reconstruct them too
    if hasattr(kmodel,'dependent_reactants'):
        new_kmodel.dependent_reactants = kmodel.dependent_reactants.copy()

    return new_kmodel
