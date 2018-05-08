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
from ..core.itemsets import make_parameter_set, make_reactant_set
from ..utils.tabdict import TabDict
from collections import namedtuple
from ..utils.namespace import *


class ReversibleMichaelisMenten(KineticMechanism):
    """A reversible uni uni reaction enzmye class"""

    Reactants = make_reactant_set(__name__, ['substrate', 'product'])

    Parameters = make_parameter_set(    __name__,
                                        {
                                        'vmax_forward':[ODE,MCA,QSSA],
                                        'k_equilibrium':[ODE,MCA,QSSA],
                                        'km_substrate':[ODE,MCA,QSSA],
                                        'km_product':[ODE,MCA,QSSA],
                                        'vmax_backward':[ODE,QSSA],
                                        'total_enzyme_concentration':[ODE,ELEMENTARY],
                                        'k1_fwd':[ODE,ELEMENTARY],
                                        'k1_bwd':[ODE,ELEMENTARY],
                                        'k2_fwd':[ODE,ELEMENTARY],
                                        'k2_bwd':[ODE,ELEMENTARY],
                                        })

    reactant_stoichiometry = {'substrate':-1,
                              'product':1}

    parameter_reactant_links = {
        'km_substrate':'substrate',
        'km_product':'product',
    }

    ElementaryReactions = namedtuple('ElementaryReactions',['r1f',
                                                            'r1b',
                                                            'r2f',
                                                            'r2b',
                                                            ])


    def __init__(self, name, reactants, parameters=None):
        # FIXME dynamic linking, separaret parametrizations from model init
        # FIXME Reaction has a mechanism, and this is a mechanism
        KineticMechanism.__init__(self, name, reactants, parameters)

    def get_qssa_rate_expression(self):
        s = self.reactants.substrate.symbol
        p = self.reactants.product.symbol

        kms = self.parameters.km_substrate.symbol
        kmp = self.parameters.km_product.symbol

        keq = self.parameters.k_equilibrium.symbol
        vmaxf = self.parameters.vmax_forward.symbol

        common_denominator = 1 + s/kms + p/kmp

        bwd_nominator = vmaxf/keq * p/kms

        fwd_nominator = vmaxf     * s/kms

        forward_rate_expression = fwd_nominator/common_denominator
        backward_rate_expression = bwd_nominator/common_denominator
        rate_expression = forward_rate_expression-backward_rate_expression

        self.reaction_rates = TabDict([('v_net', rate_expression),
                                       ('v_fwd', forward_rate_expression),
                                       ('v_bwd', backward_rate_expression),
                                       ])

        self.expressions = {s: -rate_expression,
                            p:    rate_expression}


        self.expression_parameters = self.get_parameters_from_expression(rate_expression)


    def update_qssa_rate_expression(self):
        s = self.reactants.substrate.symbol
        p = self.reactants.product.symbol

        self.expressions = {s: -self.reaction_rates['v_net'],
                            p:  self.reaction_rates['v_net']}



    def get_full_rate_expression(self):
        # Calculate rates uppon initialization
        if not hasattr(self,'rate_constants'):
            self.calculate_rate_constants()

        subs = self.reactants

        enzyme_complex = 'EC_'+self.name

        r1f = sympify(subs.substrate+"*"+self.name+'*'+'k1_fwd_'+self.name)
        r1b = sympify(enzyme_complex+'*k1_bwd_'+self.name)
        r2f = sympify(enzyme_complex+'*k2_fwd_'+self.name)
        r2b = sympify(subs.product  +"*"+self.name+'*k2_bwd_'+self.name)

        self.reaction_rates = TabDict([('r1f', r1f),
                                       ('r1b', r1b),
                                       ('r2f', r2f),
                                       ('r2b', r2b),
                                       ])

        self.expressions = {subs.substrate: r1b - r1f,
                            subs.product: r2f - r2b,
                            enzyme_complex: r1f - r1b - r2f + r2b,
                            self.name: r1b - r1f - r2b + r2f}

        self.expression_parameters = self.get_parameters_from_expression('rate_constants')




    def calculate_rate_constants(self):
        # Calcuate elementary rates

        # TODO: Check that the not None params are compatible
        # Param families/ set ?

        params = self.parameters
        k1_bwd = params.vmax_backward / params.total_enzyme_concentration
        k2_fwd = params.vmax_forward  / params.total_enzyme_concentration

        k1_fwd = (k1_bwd + k2_fwd ) / params.km_substrate
        k2_bwd = (k1_bwd + k2_fwd ) / params.km_product

        self.rate_constants = self.RateConstants( k1_fwd = k1_fwd,
                                                  k2_fwd = k2_fwd,
                                                  k1_bwd = k1_bwd,
                                                  k2_bwd = k2_bwd,
                                                  )

        # self.set_dynamic_attribute_links(self._rates)
        subs = self.reactants
        enzyme_complex = 'EC_'+self.name

        r1f = ElementrayReactionStep([self.name,subs.substrate],
                                     [enzyme_complex],
                                     'k1_fwd_'+self.name)
        r1b = ElementrayReactionStep([enzyme_complex],
                                     [self.name,subs.substrate],
                                     'k1_bwd_'+self.name)

        r2f = ElementrayReactionStep([enzyme_complex],
                                     [self.name,subs.product],
                                     'k2_fwd_'+self.name)
        r2b = ElementrayReactionStep([self.name,subs.product],
                                     [enzyme_complex],
                                     'k2_bwd_'+self.name)

        self.elementary_reactions = self.ElementaryReactions( r1f = r1f,
                                                              r1b = r1b,
                                                              r2f = r2f,
                                                              r2b = r2b,)
