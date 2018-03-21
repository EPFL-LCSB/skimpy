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
from scipy.sparse import coo_matrix
from numpy import array

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

    row = 0
    for this_variable in variables:
        column = 0
        for this_reaction in kinetic_model.reactions.values():
            substrate_names = this_reaction.mechanism.substrates._asdict().values()
            definitions = this_reaction.mechanism.substrates._asdict().keys()
            if this_variable.__str__() in substrate_names:
                N_substrate = len([1 for this_def, this_var in zip(definitions, substrate_names)
                                      if ('substrate' in this_def) and (this_variable.__str__() in this_var)])
                N_product = len([1 for this_def, this_var in zip(definitions, substrate_names)
                                      if ('product' in this_def) and (this_variable.__str__() in this_var)])
                values.append(N_substrate+N_product)
                rows.append(row)
                columns.append(column)
            column += 1
        row += 1

    shape = (len(variables), len(kinetic_model.reactions))

    stoichiometric_matrix = coo_matrix((values, (rows, columns)), shape = shape ).tocsr()

    return stoichiometric_matrix