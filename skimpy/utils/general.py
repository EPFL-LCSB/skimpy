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

from skimpy.utils.tabdict import TabDict
from skimpy.utils.namespace import *

from scipy.sparse import coo_matrix
from numpy import array
from sympy import Symbol
import re

from copy import deepcopy

def join_dicts(dicts):
    joined_dict = {}
    for dictionary in dicts:
        joined_dict.update(dictionary)

    return joined_dict


def get_stoichiometry(kinetic_model, variables):
    # get stoichometriy

    rows = []
    columns = []
    values = []

    row_ix = 0
    for this_variable in variables:
        column_ix = 0
        for this_reaction in kinetic_model.reactions.values():

            reaction_items = this_reaction.reactant_stoichiometry.items()
            reactant_names = [x.name for x,_ in reaction_items]
            if str(this_variable) in reactant_names:

                N = sum([this_stoich for this_var, this_stoich in reaction_items
                         if str(this_variable) == this_var.name])
                #Convert to real integer i
                if float(N).is_integer():
                    N = int(N)
                values.append(N)
                rows.append(row_ix)
                columns.append(column_ix)

            column_ix += 1
        row_ix += 1

    shape = (len(variables), len(kinetic_model.reactions))

    stoichiometric_matrix = coo_matrix((values, (rows, columns)), shape = shape ).tocsr()

    return stoichiometric_matrix

def check_is_symbol(s_in):
    if not isinstance(s_in, Symbol):
        return Symbol(s_in)
    else:
        return s_in


def sanitize_cobra_vars(met_name):
    # Remove dashes
    clean_met_name = re.sub(r"([a-z])\-([a-z])", r"\1_\2", str(met_name), 0, re.IGNORECASE)
    # Add underscore for variables names that start with a number
    clean_met_name = re.sub(r"(^[0-9])", r"_\1", str(clean_met_name), 0, re.IGNORECASE)
    return clean_met_name


def get_all_subclasses(cls):
    all_subclasses = []

    for subclass in cls.__subclasses__():
        all_subclasses.append(subclass)
        all_subclasses.extend(get_all_subclasses(subclass))

    return all_subclasses


def make_subclasses_dict(cls):
    the_dict = {x.__name__:x for x in get_all_subclasses(cls)}
    the_dict[cls.__name__] = cls
    return the_dict


def robust_index(in_var):
    """
    Indexing can be done with symbols or strings representing the symbol,
    so we harmonize it by returning the name of the symbol if the input is of
    type symbol

    :param in_var:
    :type in_var: str or sympy.Symbol
    :return:
    """

    if isinstance(in_var, str):
        return in_var
    elif isinstance(in_var, Symbol):
        return in_var.name
    else:
        raise TypeError('Value should be of type str or sympy.Symbol')


def get_all_reactants(model):
    reactants = TabDict([])
    for this_reaction in model.reactions.values():
        this_rectants = TabDict([(v.name,v) for v in this_reaction.mechanism.reactants.values()])

        for this_modifier in this_reaction.modifiers.values():
            this_reactants = TabDict( (r.name,r)
                                      for k,r in this_modifier.reactants.items())
            reactants.update(this_reactants)

        reactants.update(this_rectants)

    reactants = deepcopy(reactants)

    for r in reactants.values():
        r.type = VARIABLE

    return reactants