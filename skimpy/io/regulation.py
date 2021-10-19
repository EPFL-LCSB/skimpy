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
def load_enzyme_regulation(kmodel, df_regulations_all):
    """
    Test function for loading regulations into an existing kmodel

    :param kmodel: Existing kinetic model
    :param df_regulations_all: a pandas df that has the following columns: reaction_id | regulator | regulation
                                reaction_id must correspond with the names in the kinetic model
                                reactor is the name of a metabolite, with the compartment after the _
                                regulation is either 'inhibitionr', 'activation', or 'competitive inhibition'
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
        the_reactants = this_rxn.mechanism.reactants
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

            # Add existing kinetic modifiers
            modifiers = this_rxn.modifiers
            for the_modifier in modifiers:
                new_rxn.modifiers[the_modifier] = modifiers[the_modifier]

        # If there are no competitive regulations, just recreate the reaction first since the other modifications are applied on top
        else:
            TheMechanism = this_rxn.mechanism.__class__
            the_reactants = this_rxn.mechanism.reactants
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
                new_rxn.modifiers[this_regulator] = InhibitionModifier(this_regulator)
            elif this_regulation == 'activation':
                new_rxn.modifiers[this_regulator] = ActivationModifier(this_regulator)

        # Finally add the new reaction to rxns
        new_kmodel.add_reaction(new_rxn)

    # Resbuild compartments
    for the_comp in kmodel.compartments:
        new_comp = Compartment(the_comp)
        new_kmodel.add_compartment(new_comp)

    # Assing compartments
    for the_met in kmodel.reactants.values():

        met = new_kmodel.reactants[the_met.name]
        if not the_met.compartment is None:
            comp = new_kmodel.compartments[the_met.compartment]
            met.compartment = comp

    # Populate the kinmodel.parameters TabDict
    parameter_init_dict = dict()
    for rxn_obj in new_kmodel.reactions.values():
        # initalize empty param list
        parameter_init_dict[rxn_obj.name] = rxn_obj.mechanism.Parameters()

    new_kmodel.parametrize_by_reaction(parameter_init_dict)

    # Boundary Conditions
    new_kmodel.boundary_conditions = kmodel.boundary_conditions.copy()

    # # Do not forget to add the value of the BC! TODO: need to sort this out but maybe it is already sorted
    # reactant.value = the_dict['parameters'][str(reactant.symbol)]

    # Update the parameters with all pre-existing parameter values
    this_params = new_kmodel.parameters
    for parameter, symbol in kmodel.parameters.items():
        try:
            this_params[parameter].value = symbol.value
        except KeyError:
            pass

    # Initial conditions
    new_kmodel.initial_conditions = kmodel.initial_conditions.copy()

    # If the model the model has computed moieties we reconstruct them too
    if hasattr(kmodel,'dependent_reactants'):
        new_kmodel.dependent_reactants = kmodel.dependent_reactants.copy()

    return new_kmodel