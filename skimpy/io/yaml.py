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

from skimpy.utils import TabDict
from skimpy.core import Item, Reactant, Parameter, Reaction, BoundaryCondition, \
    ConstantConcentration, KineticModel, ExpressionModifier
from skimpy.mechanisms import *
from skimpy.utils.general import make_subclasses_dict
from skimpy.utils.namespace import PARAMETER, VARIABLE

ALL_MECHANISM_SUBCLASSES = make_subclasses_dict(KineticMechanism)
ALL_BOUNDARY_SUBCLASSES = make_subclasses_dict(BoundaryCondition)
ALL_MODIFIER_SUBCLASSES = make_subclasses_dict(ExpressionModifier)

#TODO We need to do better?
ALL_GENERIC_MECHANISM_SUBCLASSES = {'Convenience': make_convenience,
                                    'GeneralizedReversibleHill': make_generalized_reversible_hill_n_n,
                                    'IrrevMichaelisMenten': make_irrev_m_n_michaelis_menten,
                                    'IrrevMassaction': make_irrev_massaction,
                                    'ConvenienceInhibited': make_convenience_with_inhibition}

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
                       'initial_conditions']

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
        return dumper.represent_str(data.name)

def mechanism_representer(dumper, data):
    the_dict = {k:v.name for k,v in data.reactants.items()}
    the_dict['class'] = data.__class__.__name__
    _find = lambda s: the_dict['class'].find(s) >= 0
    if any(map(_find , ALL_GENERIC_MECHANISM_SUBCLASSES)):
        the_dict['mechanism_stoichiometry'] = data.reactant_stoichiometry

        #Clean the class name (will be reconstructed from stoichometry)
        the_dict['class'] = re_sub(data.__class__.suffix,'',the_dict['class'])

    elif any(map(_find , ALL_MODIFIER_SUBCLASSES)):
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



yaml.add_representer(TabDict, SafeRepresenter.represent_dict)
yaml.add_representer(Reactant, reactant_representer)
yaml.add_representer(Parameter, parameter_representer)
yaml.add_representer(Reaction, reaction_representer)

for the_class in ALL_MECHANISM_SUBCLASSES.values():
    yaml.add_representer(the_class, mechanism_representer)
for the_class in ALL_BOUNDARY_SUBCLASSES.values():
    yaml.add_representer(the_class, boundary_condition_representer)


def export_to_yaml(model, path=None, **kwargs):

    dict_model = vars(model)
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
    :param classdict:
    :return:
    """
    classname = classdict.pop('class')
    _find = lambda s: classname.find(s) >= 0
    try:
        if any(map(_find, ALL_MODIFIER_SUBCLASSES)):
            return ALL_MECHANISM_SUBCLASSES[classname]

        stoich_dict = classdict.pop('mechanism_stoichiometry')
        make_mechanism = get_generic_constructor(classname)
        stoichiometry = list(stoich_dict.values())
        return make_mechanism(stoichiometry)

    except KeyError:
        return ALL_MECHANISM_SUBCLASSES[classname]


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
        the_dict = yaml.load(fid)

    new = KineticModel(name = the_dict['name'])

    # Rebuild the reactions
    for the_reaction in the_dict['reactions'].values():
        TheMechanism = get_mechanism(the_reaction['mechanism'])
        the_reactants = TheMechanism.Reactants(**the_reaction['mechanism'])
        new_reaction = Reaction(name=the_reaction['name'],
                                mechanism=TheMechanism,
                                reactants=the_reactants)
        # Add kinetic modifiers
        modifiers = the_reaction['modifiers']
        for the_modifier in modifiers.values():
            TheModifier = get_mechanism(the_modifier)
            new_modifier = TheModifier(**the_modifier)
            new_reaction.modifiers[new_modifier.name] = new_modifier

        new.add_reaction(new_reaction)

    # Populate the kinmodel.parameters TabDict
    parameter_init_dict = dict()
    for rxn_obj in new.reactions.values():
        # initalize empty param list
        parameter_init_dict[rxn_obj.name] = rxn_obj.mechanism.__class__.Parameters()
    new.parametrize_by_reaction(parameter_init_dict)

    # Boundary Conditions
    for the_bc_dict in the_dict['boundary_conditions'].values():
        TheBoundaryCondition = ALL_BOUNDARY_SUBCLASSES[the_bc_dict.pop('class')]
        reactant = new.reactants[the_bc_dict.pop('reactant')]
        the_bc = TheBoundaryCondition(reactant, **the_bc_dict)
        new.add_boundary_condition(the_bc)

        # Do not forget to add the value of the BC!
        reactant.value = the_dict['parameters'][str(reactant.symbol)]


    # Parameter assignment based on what parameters have been stored
    for rxn_obj in new.reactions.values():
        # Look into parameters for assignment
        for p in rxn_obj.parameters.values():
            # Try to find the parameter in our YAML
            try:
                p.value = the_dict['parameters'][str(p.symbol)]
            except KeyError:
                #No value found
                pass

    # Initial conditions
    for the_ic, value in the_dict['initial_conditions'].items():
        new.initial_conditions[the_ic] = value

    # Do not forget to update parameters
    new.update()

    return new