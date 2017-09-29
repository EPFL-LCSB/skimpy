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

from sympy import *
from sympy.utilities.autowrap import ufuncify
from numpy import array

class ODEFunction:
    def __init__(self,variables,expr,parameters):
        self.variables  = variables
        self.expr = expr
        self.parameters  = list(parameters.values())

        #Create a binary function
        sym_vars = tuple(symbols(variables+list(parameters.keys())))
        # Awsome sympy magic
        self.function = ufuncify(sym_vars,list(expr.values()))

    def __call__(self,t,y):
        input_vars = list(y)+self.parameters
        result = self.function(*input_vars)
        return array(result)


def make_ode_fun(model,sim_type):

    # Gete all variables and expressions (Better solution with types?)
    if sim_type == 'QSSA':
        all_data = [ this_reaction.QSSA_rate_expressions() \
                        for this_reaction in model.reactions]

    elif sim_type == 'tQSSA':
        all_data = [ this_reaction.tQSSA_rate_expressions() \
                        for this_reaction in model.reactions]

    elif sim_type == 'full':
        all_data = [ this_reaction.full_rate_expressions() \
                        for this_reaction in model.reactions]

    all_vars, all_expr, all_param = list(zip(*all_data))

    # Flatten all the lists
    flatten_list = lambda this_list: [item for sublist in this_list for item in sublist]

    all_vars = flatten_list(all_vars)
    all_param = join_dicts(all_param)

    #Get unique set of all the variables
    variables = list(set(all_vars))

    expr = dict.fromkeys(variables,0.0)

    # Sum up all expressions
    for this_expr in all_expr:
        for this_expr_key in this_expr:
            expr[this_expr_key] += this_expr[this_expr_key]

    # Apply boundary conditions
    for this_boundary in model.boundaries:
        this_boundary(expr)

    # Make vector function from expressions
    return ODEFunction(variables, expr, all_param)


def join_dicts(dicts):
    joined_dict = {}
    for dictionary in dicts:
        joined_dict.update(dictionary)

    return joined_dict
