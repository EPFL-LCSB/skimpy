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

import yaml
from yaml.representer import SafeRepresenter
from re import sub as re_sub
import re

import numpy as np
from scipy.sparse import csr_matrix, csc_matrix

from skimpy.utils import TabDict, iterable_to_tabdict
from skimpy.core import Item, Reactant, Parameter, Reaction, BoundaryCondition, \
    ConstantConcentration, KineticModel, ExpressionModifier
from skimpy.core.compartments import Compartment
from skimpy.mechanisms import *
from skimpy.utils.general import make_subclasses_dict, get_stoichiometry, get_all_reactants
from skimpy.utils.namespace import PARAMETER, VARIABLE


from skimpy.analysis.mca.utils import get_reduced_stoichiometry

def get_mechanism_subclasses():
    return make_subclasses_dict(KineticMechanism)

def get_boundary_subclasses():
    return make_subclasses_dict(BoundaryCondition)

def get_modifier_subclasses():
    return make_subclasses_dict(ExpressionModifier)

#TODO We need to do better?
ALL_GENERIC_MECHANISM_SUBCLASSES = TabDict([
('Convenience', make_convenience),
('GeneralizedReversibleHill', make_generalized_reversible_hill_n_n),
('H1GeneralizedReversibleHill',make_generalized_reversible_hill_n_n_h1),
('IrrevMichaelisMenten', make_irrev_m_n_michaelis_menten),
('IrrevHillNM', make_irrev_m_n_hill),
('IrrevMassaction', make_irrev_massaction),
('RevMassaction', make_rev_massaction),
('ConvenienceInhibited', make_convenience_with_inhibition),
('H1GeneralizedReversibleHillInhibited', make_generalized_reversible_hill_n_n_h1_with_inhibition),
])

FIELDS_TO_SERIALIZE = [
                       # 'variables',
                       # 'ode_fun',
                       'boundary_conditions',
                       'parameters',
                       #'logger',
                       # '_simtype',
                       'constraints',
                       # 'solver',
                       # '_recompiled', '_modifed',
                       'reactions',
                       # '_modified',
                       'name',
                       'initial_conditions',
                       'dependent_reactants',
                       'compartments',
                       'reactants',
                        ]

#----------------------------------------------------------------
#                       Model serialization
#----------------------------------------------------------------


def parameter_representer(dumper, data):
    if data.value is not None:
        return dumper.represent_float(data.value)
    else:
        return  dumper.represent_none(data.value)

def reactant_representer(dumper, data):
    if data.type == PARAMETER:
        # Store it as a parameter
        return parameter_representer(dumper,data)
    elif data.type == VARIABLE:
        dict = {'name':data.name,
                'compartment': None if data.compartment is None else data.compartment.name}
        return dumper.represent_dict(dict)

def compartment_representer(dumper, data):
    return dumper.represent_str(data.name)

def mechanism_representer(dumper, data):
    the_dict = {k:v.name for k,v in data.reactants.items()}
    the_dict['class'] = data.__class__.__name__

    _find = lambda s: the_dict['class'].find(s) >= 0
    if any(map(_find , ALL_GENERIC_MECHANISM_SUBCLASSES)):
        the_dict['mechanism_stoichiometry'] = data.reactant_stoichiometry

        #Clean the class name (will be reconstructed from stoichometry)
        the_dict['class'] = re_sub(data.__class__.suffix,'',the_dict['class'])

    elif any(map(_find , get_modifier_subclasses())):
        the_dict['mechanism_stoichiometry'] = data.reactant_stoichiometry

    return dumper.represent_dict(the_dict)

def reaction_representer(dumper,data):
    return dumper.represent_dict(vars(data))

def boundary_condition_representer(dumper, data):
    the_dict = dict()
    classname = data.__class__.__name__
    the_dict['class'] = classname
    the_dict['reactant'] = data.reactant.name
    if classname == ConstantConcentration.__name__:
        # There is no extra fields
        pass
    else:
        raise NotImplementedError('This modifier class is not implemented yet')
        the_dict['value'] = the_value

    return dumper.represent_dict(the_dict)


def refresh_representers():
    """
    This function refreshes the representers of custom-made reaction mechanisms 
    (n-to-m Convenience Kinetics for example)
    """
    yaml.add_representer(TabDict, SafeRepresenter.represent_dict)
    yaml.add_representer(Reactant, reactant_representer)
    yaml.add_representer(Parameter, parameter_representer)
    yaml.add_representer(Reaction, reaction_representer)
    yaml.add_representer(Compartment, compartment_representer)

    for the_class in get_mechanism_subclasses().values():
        yaml.add_representer(the_class, mechanism_representer)
    for the_class in get_boundary_subclasses().values():
        yaml.add_representer(the_class, boundary_condition_representer)


def export_to_yaml(model, path=None, **kwargs):

    # In case new custom mechanisms have been made 
    refresh_representers()

    dict_model = vars(model)
    # Add parameters that are properties
    dict_model['parameters'] = model.parameters
    dict_model['reactants'] = get_all_reactants(model)
    fields_not_to_serialize = [x for x in dict_model if not x in FIELDS_TO_SERIALIZE]
    [dict_model.pop(k) for k in fields_not_to_serialize]

    if path is not None:
        with open(path, 'w') as fid:
            yaml.dump(dict_model, fid, default_flow_style=False)
            return True
    else:
        return yaml.dump(dict_model, default_flow_style=False)


#----------------------------------------------------------------
#                       Model loading
#----------------------------------------------------------------

def get_mechanism(classdict):
    #TODO Make nice and more readable
    """
    This function should construct mechanism for the generic types
    and get the mechanism for modifiers and the defined types

    :param classdict: Looks like {'class':'MichaelisMenten',
                                  'substrate1':'atp_c',
                                  'product1'  :'gtp_c',
                                  }

    :return:
    """
    classname = classdict.pop('class')

    # Checks if the classname is a modifier (Modifiers are subclasses of Mechanisms)
    # SPLIT FOR MECHANISM vs Modifier
    _find = lambda s: classname.find(s) >= 0
    try:
        if any(map(_find, get_modifier_subclasses())):
            # If the class name is indeed a modifier, get the actual
            return get_mechanism_subclasses()[classname]

        # TODO Make this realiable and nice !!!!
        # TODO the way we create the mechanism can be easily
        stoich_dict = classdict.pop('mechanism_stoichiometry')
        make_mechanism = get_generic_constructor(classname)

        index_stoich = [(int(re.findall(r'\d+',k)[0]), v)
                        for k, v in stoich_dict.items()]

        stoichiometry = [v for k,v in sorted(index_stoich)]
        return make_mechanism(stoichiometry)

    except KeyError:
        return get_mechanism_subclasses()[classname]


def get_generic_constructor(s):
    for name, constructor in ALL_GENERIC_MECHANISM_SUBCLASSES.items():
        if s.startswith(name):
            return constructor

def get_stoich(s):
    #TODO can we generalize this so it can include inhibitors?
    #Using regexp on something like e.g. _s1_s1_p1_i1_
    splitted = s.split('_')
    return [int(x.replace('m','-')) for x in splitted[1:]]

def load_yaml_model(path):
    with open(path,'r') as fid:
        the_dict = yaml.full_load(fid)

    new = KineticModel(name = the_dict['name'])

    # Rebuild the reactions
    for the_reaction in the_dict['reactions'].values():
        TheMechanism = get_mechanism(the_reaction['mechanism'])
        the_reactants = TheMechanism.Reactants(**the_reaction['mechanism'])
        the_enzyme = the_reaction['enzyme']

        new_reaction = Reaction(name=the_reaction['name'],
                                mechanism=TheMechanism,
                                reactants=the_reactants,
                                enzyme = the_enzyme,
                                )
        # Add kinetic modifiers
        modifiers = the_reaction['modifiers']
        for the_mod_name, the_modifier in modifiers.items():
            TheModifier = get_mechanism(the_modifier)
            new_modifier = TheModifier(**the_modifier, name=the_mod_name, reaction=new_reaction)
            new_reaction.modifiers[new_modifier.name] = new_modifier

        new.add_reaction(new_reaction)

    #Resbuild compartments
    for the_comp in the_dict['compartments'].values():
        new_comp = Compartment(the_comp)
        new.add_compartment(new_comp)

    # Assing compartments
    for the_met in the_dict['reactants'].values():
        met = new.reactants[the_met['name']]
        if not the_met['compartment'] is None:
            comp = new.compartments[the_met['compartment']]
            met.compartment = comp

    # Populate the kinmodel.parameters TabDict
    parameter_init_dict = dict()
    for rxn_obj in new.reactions.values():
        # initalize empty param list
        parameter_init_dict[rxn_obj.name] = rxn_obj.mechanism.__class__.Parameters()

    new.parametrize_by_reaction(parameter_init_dict)

    # Boundary Conditions
    for the_bc_dict in the_dict['boundary_conditions'].values():
        TheBoundaryCondition = get_boundary_subclasses()[the_bc_dict.pop('class')]
        reactant = new.reactants[the_bc_dict.pop('reactant')]
        the_bc = TheBoundaryCondition(reactant, **the_bc_dict)
        new.add_boundary_condition(the_bc)

        # Do not forget to add the value of the BC!
        reactant.value = the_dict['parameters'][str(reactant.symbol)]


    # Parameter assignment based on what parameters have been stored
    # for rxn_obj in new.reactions.values():
    #     # Look into parameters for assignment
    #     for p in rxn_obj.parameters.values():
    #         # Try to find the parameter in our YAML
    #         try:
    #             p.value = the_dict['parameters'][str(p.symbol)]
    #         except KeyError:
    #             #No value found
    #             pass
    this_params = new.parameters
    for parameter, value in the_dict['parameters'].items():
        try:
            this_params[parameter].value = value
        except KeyError:
            pass

    # Initial conditions
    for the_ic, value in the_dict['initial_conditions'].items():
        new.initial_conditions[the_ic] = value


    #If the model the model has computed moieties we reconstruct them too
    if 'dependent_reactants' in the_dict.keys():
        rebuild_dependent_mets(new, the_dict)

    return new


def rebuild_dependent_mets(new, the_dict):
    new.dependent_reactants = iterable_to_tabdict([new.reactants[r]
                                                   for r in the_dict["dependent_reactants"]])
    new.dependent_variables_ix = [i for i, e in enumerate(new.reactants.keys())
                                  if e in new.dependent_reactants.keys()]
    new.independent_variables_ix = [i for i, e in enumerate(new.reactants.keys())
                                    if e not in new.dependent_reactants.keys()]

    new.variables = TabDict([(k, v.symbol) for k, v in new.reactants.items()])

    reduced_stoichiometry, conservation_relation, _, _ \
        = get_reduced_stoichiometry(new,
                                    new.variables,
                                    all_dependent_ix=new.dependent_variables_ix)

    new.reduced_stoichiometry = reduced_stoichiometry
    new.conservation_relation = conservation_relation
