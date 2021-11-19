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
from numpy import abs as np_abs

from .mechanism import KineticMechanism,ElementrayReactionStep
from ..core.reactions import Reaction
from ..utils.tabdict import TabDict
from collections import namedtuple
from ..core.itemsets import make_parameter_set, make_reactant_set
from skimpy.utils.general import make_subclasses_dict
from ..utils.namespace import *
from .utils import stringify_stoichiometry


def make_generalized_reversible_hill_n_n(stoichiometry):

    """

    :param stoichiometry is a list of the reaction stoichioemtry
    """

    # This refresh all subclasses and fetches already create mechanism classes
    ALL_MECHANISM_SUBCLASSES = make_subclasses_dict(KineticMechanism)

    new_class_name = "GeneralizedReversibleHill"\
                     + "_{0}".format(stringify_stoichiometry(stoichiometry))

    if new_class_name in ALL_MECHANISM_SUBCLASSES.keys():
        return ALL_MECHANISM_SUBCLASSES[new_class_name]

    class GeneralizedReversibleHill(KineticMechanism):
        """
        A reversible hill N-N enzyme class
        e.g.: A + B -> C + D

        """

        # This will allow to abuse this kinetic for practical purposes
        # #check_if stoichometry_is one
        # if any(np_abs(stoichiometry) > 1)  \
        #    or not (stoichiometry.count(-1) == stoichiometry.count(1)):
        #     raise ValueError('Stoichiometry needs to be 1 and n to n substrates! ')

        suffix = "_{0}".format(stringify_stoichiometry(stoichiometry))

        reactant_list = []
        parameter_list = {'vmax_forward': [ODE, MCA, QSSA],
                          'kcat_forward':[ODE,MCA,QSSA],
                          'k_equilibrium': [ODE, MCA, QSSA],
                          'hill_coefficient': [ODE, MCA, QSSA]}

        parameter_reactant_links = {}
        reactant_parameter_links = {}
        reactant_stoichiometry = {}

        num_substrates = 1
        num_products = 1
        for s in stoichiometry:
            if s < 0:
                substrate = 'substrate{}'.format(num_substrates)
                km_substrate ='km_substrate{}'.format(num_substrates)
                reactant_list.append(substrate)
                parameter_list[km_substrate] = [ODE, MCA, QSSA]
                parameter_reactant_links[km_substrate] = substrate
                reactant_parameter_links[substrate] = km_substrate
                reactant_stoichiometry[substrate] = float(s)
                num_substrates += 1

            if s > 0:
                product = 'product{}'.format(num_products)
                km_product ='km_product{}'.format(num_products)
                reactant_list.append(product)
                parameter_list[km_product] = [ODE, MCA, QSSA]
                parameter_reactant_links[km_product] = product
                reactant_parameter_links[product] = km_product
                reactant_stoichiometry[product] = float(s)
                num_products += 1

        Reactants = make_reactant_set(__name__ + suffix, reactant_list)

        Parameters = make_parameter_set(__name__ + suffix, parameter_list)

        ElementaryReactions = namedtuple('ElementaryReactions',[])


        def __init__(self, name, reactants, parameters=None, **kwargs):
            KineticMechanism.__init__(self, name, reactants, parameters, **kwargs)

        def get_qssa_rate_expression(self):
            # NOTE: THis should not be done based on the symbol of the reactant
            # In case a reactant is duplicate due to mechanistic reasons!
            # reactant_km_relation = {self.reactants[v].symbol: k
            #                         for k, v in self.parameter_reactant_links.items()}

            substrates = TabDict([(k, self.reactants[k])
                                  for k in self.reactant_list
                                  if k.startswith('substrate')])

            products = TabDict([(k, self.reactants[k])
                                for k in self.reactant_list
                                if k.startswith('product')])


            keq = self.parameters.k_equilibrium.symbol
            if self.enzyme is None:
                vmaxf = self.parameters.vmax_forward.symbol
            else:
                vmaxf = self.parameters.kcat_forward.symbol * \
                        self.reactants.enzyme.symbol

            fwd_nominator = vmaxf
            bwd_nominator = vmaxf/keq

            for type, this_substrate in substrates.items():
                s = this_substrate.symbol
                kms = self.parameters[self.reactant_parameter_links[type]].symbol
                stoich = self.reactant_stoichiometry[type]
                fwd_nominator *= (s/kms)**abs(stoich)
                bwd_nominator *= kms**(-1*abs(stoich))

            for type, this_product in products.items():
                p = this_product.symbol
                stoich = self.reactant_stoichiometry[type]
                bwd_nominator *= p**abs(stoich)

            common_denominator = 1

            h = self.parameters.hill_coefficient.symbol
            for this_product,this_substrate in zip(substrates.items(),products.items()):
                substrate_type = this_substrate[0]
                product_type = this_product[0]
                s = this_substrate[1].symbol
                p = this_product[1].symbol
                kms = self.parameters[self.reactant_parameter_links[substrate_type]].symbol
                kmp = self.parameters[self.reactant_parameter_links[product_type]].symbol

                common_denominator *= (1 + (p/kmp + s/kms)**h)\
                                      / (p/kmp + s/kms)**(h-1)

            forward_rate_expression = fwd_nominator/common_denominator
            backward_rate_expression = bwd_nominator/common_denominator
            rate_expression = forward_rate_expression-backward_rate_expression

            self.reaction_rates = TabDict([('v_net', rate_expression),
                                           ('v_fwd', forward_rate_expression),
                                           ('v_bwd', backward_rate_expression),
                                           ])

            expressions = {}

            # TODO Find a better solution to handle duplicate substrates
            # The dict currently does not allow for this 
            for type, this_substrate in substrates.items():
                s = this_substrate.symbol
                stoich = self.reactant_stoichiometry[type]
                if s in expressions.keys():
                    expressions[s] += stoich * rate_expression
                else:
                    expressions[s] = stoich*rate_expression

            for type, this_product in products.items():
                p = this_product.symbol
                stoich = self.reactant_stoichiometry[type]
                if p in expressions.keys():
                    expressions[p] += stoich * rate_expression
                else:
                    expressions[p] = stoich*rate_expression

            self.expressions = expressions
            self.expression_parameters = self.get_parameters_from_expression(rate_expression)

        def update_qssa_rate_expression(self):

            substrates = {k:r for k,r in self.reactants.items()
                          if k.startswith('substrate')}

            products= {k:r for k,r in self.reactants.items()
                          if k.startswith('product')}

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

    GeneralizedReversibleHill.__name__ += GeneralizedReversibleHill.suffix

    return GeneralizedReversibleHill