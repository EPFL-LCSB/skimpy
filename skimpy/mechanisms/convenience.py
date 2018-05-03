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


from sympy import sympify
from .mechanism import KineticMechanism,ElementrayReactionStep
from ..core.reactions import Reaction
from ..utils.tabdict import TabDict
from collections import namedtuple
from ..core.itemsets import make_parameter_set, make_reactant_set
from ..utils.namespace import *


def make_convenience(stoichioemtry):

    """

    :param stoichioemtry is a list of the reaction stoichioemtry
    """
    class Convenience(KineticMechanism):
        """A reversible N-M enyme class """

        suffix = "_{0}".format(stoichioemtry)

        reactant_list = []
        parameter_list = {'vmax_forward': [ODE, MCA, QSSA],
                          'k_equilibrium': [ODE, MCA, QSSA], }

        parameter_reactant_links = {}
        reactant_stoichiometry = {}

        num_substrates = 1
        num_products = 1
        for s in stoichioemtry:
            if s < 0:
                substrate = 'substrate{}'.format(num_substrates)
                km_substrate ='km_substrate{}'.format(num_substrates)
                reactant_list.append(substrate)
                parameter_list[km_substrate] = [ODE, MCA, QSSA]
                parameter_reactant_links[km_substrate] = substrate
                reactant_stoichiometry[substrate] = s
                num_substrates += 1

            if s > 0:
                product = 'product{}'.format(num_products)
                km_product ='km_product{}'.format(num_products)
                reactant_list.append(product)
                parameter_list[km_product] = [ODE, MCA, QSSA]
                parameter_reactant_links[km_product] = product
                reactant_stoichiometry[product] = s
                num_products += 1

        Reactants = make_reactant_set(__name__ + suffix, reactant_list)

        Parameters = make_parameter_set(__name__ + suffix, parameter_list)

        ElementaryReactions = namedtuple('ElementaryReactions',[])


        def __init__(self, name, reactants, parameters=None):
            # FIXME dynamic linking, separaret parametrizations from model init
            # FIXME Reaction has a mechanism, and this is a mechanism
            KineticMechanism.__init__(self, name, reactants, parameters)

        def get_qssa_rate_expression(self):

            substrates = {k:r for k,r in self.reactants.items()
                          if k.startswith('substrate')}

            products= {k:r for k,r in self.reactants.items()
                          if k.startswith('product')}


            keq = self.parameters.k_equilibrium.symbol
            vmaxf = self.parameters.vmax_forward.symbol

            common_denominator_substrates = 1
            fwd_nominator = vmaxf
            bwd_nominator = vmaxf/keq

            for type, this_substrate in substrates.items():
                common_denominator_this_substrate = 1
                s = this_substrate.symbol
                kms = this_substrate.hook.symbol
                stoich = self.reactant_stoichiometry[type]
                for alpha in range(abs(stoich)):
                    common_denominator_this_substrate += (s/kms)**alpha
                # Multiply for every substrate
                common_denominator_substrates *= common_denominator_this_substrate

                fwd_nominator *= (s/kms)**abs(stoich)
                bwd_nominator *= kms**(-1*abs(stoich))

            common_denominator_products = 1
            for type, this_product in products.items():
                common_denominator_this_product = 1
                p = this_product.symbol
                kmp = this_product.hook.symbol
                stoich = self.reactant_stoichiometry[type]
                for beta in range(abs(stoich)):
                    common_denominator_this_product += (p/kmp)**beta
                # Multiply for every product
                common_denominator_products *= common_denominator_this_product

                bwd_nominator *= p**abs(stoich)


            common_denominator = common_denominator_substrates +\
                                 common_denominator_products - 1

            forward_rate_expression = fwd_nominator/common_denominator
            backward_rate_expression = bwd_nominator/common_denominator
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

            parameters = self.get_parameters_from_expression(rate_expression)
            return expressions, parameters

        """"
        Convenience kinetics has no detailed mechanism 
        """
        def get_full_rate_expression(self):
            raise NotImplementedError

        def calculate_rate_constants(self):
            raise NotImplementedError

    Convenience.__name__ += Convenience.suffix

    return Convenience