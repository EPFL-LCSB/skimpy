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


class RandBiBiReversibleMichaelisMenten(KineticMechanism):
    """A reversible random ordered bi bi reaction enzmye class"""

    Reactants = make_reactant_set(__name__, [ 'substrate1',
                                              'substrate2',
                                              'product1',
                                              'product2'])

    Parameters = make_parameter_set(__name__,
                                    {   'vmax_forward':[ODE,MCA,QSSA],
                                        'kcat_forward':[ODE,MCA,QSSA],
                                        'k_equilibrium':[ODE,MCA,QSSA],
                                        'ki_substrate1':[ODE,MCA,QSSA],
                                        'ki_substrate2':[ODE,MCA,QSSA],
                                        'km_substrate2':[ODE,MCA,QSSA],
                                        'km_product1':[ODE,MCA,QSSA],
                                        'ki_product1':[ODE,MCA,QSSA],
                                        'ki_product2':[ODE,MCA,QSSA],
                                        'vmax_backward':[ODE,QSSA],
                                        'total_enzyme_concentration':[ODE,ELEMENTARY],
                                    })

    reactant_stoichiometry= {'substrate1':-1,
                             'substrate2':-1,
                             'product1':1,
                             'product2':1}

    parameter_reactant_links = {
        'ki_substrate1':'substrate1',
        'ki_substrate2':'substrate2',
        'km_substrate2':'substrate2',
        'km_product1':'product1',
        'ki_product1':'product1',
        'ki_product2':'product2',
    }



    ElementaryReactions = namedtuple('ElementaryReactions',[])


    def __init__(self, name, reactants, parameters=None, **kwargs):
        # FIXME dynamic linking, separaret parametrizations from model init
        # FIXME Reaction has a mechanism, and this is a mechanism
        KineticMechanism.__init__(self, name, reactants, parameters, **kwargs)


    def get_qssa_rate_expression(self):

        s1 = self.reactants.substrate1.symbol
        s2 = self.reactants.substrate2.symbol
        p1 = self.reactants.product1.symbol
        p2 = self.reactants.product2.symbol

        kis1 = self.parameters.ki_substrate1.symbol
        kis2 = self.parameters.ki_substrate2.symbol
        kms2 = self.parameters.km_substrate2.symbol

        kmp1 = self.parameters.km_product1.symbol
        kip1 = self.parameters.ki_product1.symbol
        kip2 = self.parameters.ki_product2.symbol

        keq = self.parameters.k_equilibrium.symbol

        if self.enzyme is None:
            vmaxf = self.parameters.vmax_forward.symbol
        else:
            vmaxf = self.parameters.kcat_forward.symbol * \
                    self.reactants.enzyme.symbol
            
        common_denominator = 1 + s1/kis1 + s2/kis2 + p1/kip1 + p2/kip2 + \
                            s1*s2/(kis1*kms2) + p1*p2/(kip2*kmp1)

        bwd_nominator = vmaxf/keq * (p1*p2)/(kis1*kms2)

        fwd_nominator = vmaxf * (s1*s2)/(kis1*kms2)

        forward_rate_expression = fwd_nominator/common_denominator
        backward_rate_expression = bwd_nominator/common_denominator
        rate_expression = forward_rate_expression-backward_rate_expression

        self.reaction_rates = TabDict([('v_net', rate_expression),
                                       ('v_fwd', forward_rate_expression),
                                       ('v_bwd', backward_rate_expression),
                                       ])

        self.expressions = {s1: -rate_expression,
                            s2: -rate_expression,
                            p1:    rate_expression,
                            p2:    rate_expression
                            }

        self.expression_parameters = self.get_parameters_from_expression(rate_expression)

    def update_qssa_rate_expression(self):
        s1 = self.reactants.substrate1.symbol
        s2 = self.reactants.substrate2.symbol
        p1 = self.reactants.product1.symbol
        p2 = self.reactants.product2.symbol

        self.expressions = {s1: -self.reaction_rates['v_net'],
                            s2: -self.reaction_rates['v_net'],
                            p1: self.reaction_rates['v_net'],
                            p2: self.reaction_rates['v_net']}

    def get_full_rate_expression(self):
        raise NotImplementedError

    def calculate_rate_constants(self):
        raise NotImplementedError
