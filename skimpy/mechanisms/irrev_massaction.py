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

import sympy

from sympy import sympify
from .mechanism import KineticMechanism,ElementrayReactionStep
from ..core.reactions import Reaction
from ..utils.tabdict import TabDict
from collections import namedtuple
from ..core.itemsets import make_parameter_set, make_reactant_set
from skimpy.utils.general import make_subclasses_dict
from ..utils.namespace import *
from .utils import stringify_stoichiometry



def make_irrev_massaction(stoichiometry):

    """

    :param stoichiometry is a list of the reaction stoichioemtry
    """
    # This refresh all subclasses and fetches already create mechanism classes
    ALL_MECHANISM_SUBCLASSES = make_subclasses_dict(KineticMechanism)

    new_class_name = "IrrevMassaction"\
                     + "_{0}".format(stringify_stoichiometry(stoichiometry))

    if new_class_name in ALL_MECHANISM_SUBCLASSES.keys():
        return ALL_MECHANISM_SUBCLASSES[new_class_name]

    class IrrevMassaction(KineticMechanism):
        """A reversible N-M enyme class """

        suffix = "_{0}".format(stringify_stoichiometry(stoichiometry))

        reactant_list = []
        parameter_list = {'vmax_forward': [ODE, MCA, QSSA],}

        parameter_reactant_links = {}
        reactant_stoichiometry = {}

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

        ElementaryReactions = namedtuple('ElementaryReactions',[])


        def __init__(self, name, reactants, parameters=None, **kwargs):
            # FIXME dynamic linking, separaret parametrizations from model init
            # FIXME Reaction has a mechanism, and this is a mechanism
            KineticMechanism.__init__(self, name, reactants, parameters, **kwargs)

        def get_qssa_rate_expression(self):

            substrates = {k:r for k,r in self.reactants.items()
                          if k.startswith('substrate')}

            products= {k:r for k,r in self.reactants.items()
                          if k.startswith('product')}

            kf = self.parameters.vmax_forward.symbol

            forward_rate_expression = kf
            backward_rate_expression = sympy.S.Zero

            for type, this_substrate in substrates.items():
                s = this_substrate.symbol
                forward_rate_expression *= s

            rate_expression = forward_rate_expression-backward_rate_expression

            self.reaction_rates = TabDict([('v_net', rate_expression),
                                           ('v_fwd', forward_rate_expression),
                                           ('v_bwd', backward_rate_expression),
                                           ])

            expressions = {}

            for type, this_substrate in substrates.items():
                s = this_substrate.symbol
                stoich = self.reactant_stoichiometry[type]
                expressions[s] = stoich*rate_expression

            for type, this_product in products.items():
                p = this_product.symbol
                stoich = self.reactant_stoichiometry[type]
                expressions[p] = stoich * rate_expression

            self.expressions = expressions
            self.expression_parameters = self.get_parameters_from_expression(rate_expression)

        def update_qssa_rate_expression(self):

            substrates = {k:r for k,r in self.reactants.items()
                          if k.startswith('substrate')}

            products= {k:r for k,r in self.reactants.items()
                          if k.startswith('product')}
            
            for type, this_substrate in substrates.items():
                s = this_substrate.symbol
                stoich = self.reactant_stoichiometry[type]
                self.expressions[s] = stoich*self.reaction_rates['v_net']

            for type, this_product in products.items():
                p = this_product.symbol
                stoich = self.reactant_stoichiometry[type]
                self.expressions[p] = stoich*self.reaction_rates['v_net']


        """"
        Not implemented yet
        """
        def get_full_rate_expression(self):
            raise NotImplementedError

        def calculate_rate_constants(self):
            raise NotImplementedError

    IrrevMassaction.__name__ += IrrevMassaction.suffix

    return IrrevMassaction