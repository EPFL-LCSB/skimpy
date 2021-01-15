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


class UniBiReversibleHill(KineticMechanism):
    """
    A reversible hill Uni-Bi enzyme class
    e.g.: A -> B + C

    """

    reactant_list = ['substrate',
                     'product1',
                     'product2',
                     ]

    Reactants = make_reactant_set(__name__, reactant_list)

    Parameters = make_parameter_set(__name__,
                                        {
                                            'vmax_forward': [ODE, MCA, QSSA],
                                            'kcat_forward':[ODE,MCA,QSSA],
                                            'k_equilibrium':[ODE, MCA, QSSA],
                                            'hill_coefficient': [ODE, MCA, QSSA],
                                            'km_substrate': [ODE, MCA, QSSA],
                                            'km_product1':[ODE,MCA,QSSA],
                                            'km_product2': [ODE, MCA, QSSA],
                                            'vmax_backward':[ODE,QSSA],
                                        })


    reactant_stoichiometry = {'substrate': -1,
                              'product1': 1,
                              'product2': 1
                              }

    parameter_reactant_links = {
        'km_substrate':    'substrate',
        'km_product1':     'product1',
        'km_product2':     'product2',
    }


    ElementaryReactions = namedtuple('ElementaryReactions',[])


    def __init__(self, name, reactants, parameters=None, **kwargs):
        KineticMechanism.__init__(self, name, reactants, parameters, **kwargs)

    def get_qssa_rate_expression(self):
        reactant_km_relation = {self.reactants[v].symbol: k
                                for k, v in self.parameter_reactant_links.items()}

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

        s = substrates['substrate'].symbol
        p1 = products['product1'].symbol
        p2 = products['product2'].symbol

        kms = self.parameters[reactant_km_relation[s]].symbol
        kmp1 = self.parameters[reactant_km_relation[p1]].symbol
        kmp2 = self.parameters[reactant_km_relation[p2]].symbol

        h = self.parameters.hill_coefficient.symbol

        hill_effect = (s/kms+p1/kmp1*p2/kmp2)**(h-1)

        fwd_nominator = vmaxf*s/kms*hill_effect
        bwd_nominator = vmaxf/(keq*kms)*p1*p2*hill_effect

        common_denominator = 1 + (s/kms+p1/kmp1)**h \
                             + (s/kms+p2/kmp2)**h \
                             + (s/kms+p1/kmp1*p2/kmp2)**h \
                             - 2*s/kms**h

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
    GeneralizedReversibleHill kinetics has no detailed mechanism 
    """
    def get_full_rate_expression(self):
        raise NotImplementedError

    def calculate_rate_constants(self):
        raise NotImplementedError
