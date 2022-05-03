# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2022 Laboratory of Computational Systems Biotechnology (LCSB),
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

from sympy import sympify
from numpy import abs as np_abs

from .mechanism import KineticMechanism,ElementrayReactionStep
from ..core.reactions import Reaction
from ..utils.tabdict import TabDict
from collections import namedtuple
from ..core.itemsets import make_parameter_set, make_reactant_set
from skimpy.utils.general import make_subclasses_dict
from ..utils.namespace import *
from .utils import stringify_stoichiometry

import uuid
from sympy.parsing.sympy_parser import parse_expr
from sympy.parsing.sympy_parser import standard_transformations,\
    implicit_multiplication_application

TRANSFORMATIONS = (standard_transformations +(implicit_multiplication_application,))


def make_mechanism_from_expression(stoichiometry, reactants, rate_expression):
    """
    stoichiometry: integer array [-1, 1]
    parameters: list of strings ['k1', 'Keq']
    rate_expression: String expression 'k1*A*B'
    """

    # This refresh all subclasses and fetches already create mechanism classes
    ALL_MECHANISM_SUBCLASSES = make_subclasses_dict(KineticMechanism)

    # Generate a unique class name
    new_class_name = "ExpresssionBased"\
                     + "_{0}_{1}".format(stringify_stoichiometry(stoichiometry), uuid.uuid4() )

    if new_class_name in ALL_MECHANISM_SUBCLASSES.keys():
        return ALL_MECHANISM_SUBCLASSES[new_class_name]

    class ExpressionBasedKinetics(KineticMechanism):
        """
        A class building abitrary kinetics based on the expression
        """

        suffix = "_{0}".format(stringify_stoichiometry(stoichiometry))

        dictionary = {r: '1' for r in reactants}
        temp = rate_expression
        for key in dictionary.keys():
            temp = temp.replace(key, dictionary[key])

        temp_expr = parse_expression(temp)
        parameters = [str(s) for s in temp_expr.free_symbols]

        reactant_list = []
        parameter_list = {p: [ODE, ] for p in parameters}
        reactant_stoichiometry = {}

        parameter_reactant_links = {}
        reactant_parameter_links = {}

        num_substrates = 1
        num_products = 1

        for s in stoichiometry:
            if s < 0:
                substrate = 'substrate{}'.format(num_substrates)
                reactant_list.append(substrate)
                reactant_stoichiometry[substrate] = float(s)
                num_substrates += 1

            if s > 0:
                product = 'product{}'.format(num_products)
                reactant_list.append(product)
                reactant_stoichiometry[product] = float(s)
                num_products += 1

        Reactants = make_reactant_set(__name__ + suffix, reactant_list)

        Parameters = make_parameter_set(__name__ + suffix, parameter_list)

        ElementaryReactions = namedtuple('ElementaryReactions', [])

        def __init__(self, name, reactants, parameters=None, **kwargs):
            KineticMechanism.__init__(self, name, reactants, parameters, **kwargs)

        def get_qssa_rate_expression(self):
            substrates = {k: self.reactants[k] for k, s in self.reactant_stoichiometry.items()
                          if s < 0}

            products = {k: self.reactants[k] for k, s in self.reactant_stoichiometry.items()
                        if s > 0}

            parsed_rate_expression = parse_expression(rate_expression)

            self.reaction_rates = TabDict([('v_net', parsed_rate_expression),
                                           ('v_fwd', parsed_rate_expression),
                                           ('v_bwd', 0),
                                           ])

            expressions = {}

            # TODO Find a better solution to handle duplicate substrates
            # The dict currently does not allow for this
            for type, this_substrate in substrates.items():
                s = this_substrate.symbol
                stoich = self.reactant_stoichiometry[type]
                if s in expressions.keys():
                    expressions[s] += stoich * parsed_rate_expression
                else:
                    expressions[s] = stoich*parsed_rate_expression

            for type, this_product in products.items():
                p = this_product.symbol
                stoich = self.reactant_stoichiometry[type]
                if p in expressions.keys():
                    expressions[p] += stoich * parsed_rate_expression
                else:
                    expressions[p] = stoich*parsed_rate_expression

            self.expressions = expressions
            self.expression_parameters = self.get_parameters_from_expression(parsed_rate_expression)

        def update_qssa_rate_expression(self):

            substrates = {k: self.reactants[k] for k,s in self.reactant_stoichiometry.items()
                          if s < 0}

            products= {k: self.reactants[k] for k,s in self.reactant_stoichiometry.items()
                          if s > 0 }

            expressions = {}
            for type, this_substrate in substrates.items():
                s = this_substrate.symbol
                stoich = self.reactant_stoichiometry[type]
                if s in expressions.keys():
                    expressions[s] += stoich * self.reaction_rates['v_net']
                else:
                    expressions[s] = stoich*self.reaction_rates['v_net']

            for type, this_product in products.items():
                p = this_product.symbol
                stoich = self.reactant_stoichiometry[type]
                if p in expressions.keys():
                    expressions[p] += stoich * self.reaction_rates['v_net']
                else:
                    expressions[p] = stoich * self.reaction_rates['v_net']

                self.expressions = expressions


        """"
        GeneralizedReversibleHill kinetics has no detailed mechanism 
        """
        def get_full_rate_expression(self):
            raise NotImplementedError

        def calculate_rate_constants(self):
            raise NotImplementedError

    ExpressionBasedKinetics.__name__ += ExpressionBasedKinetics.suffix

    return ExpressionBasedKinetics


def parse_expression(string_expression):
    """
    Parse a string into an expression
    This is the reason this class is very fragile
    """
    exp = parse_expr(string_expression, transformations=TRANSFORMATIONS)

    return exp



