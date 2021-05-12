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
from ..core.itemsets import make_parameter_set, make_reactant_set, Reactant
from ..utils.tabdict import TabDict
from collections import namedtuple
from ..utils.namespace import *


class ReversibleMichaelisMenten(KineticMechanism):
    """A reversible uni uni reaction enzmye class"""

    Reactants = make_reactant_set(__name__, ['substrate', 'product'])

    Parameters = make_parameter_set(__name__,
                                    {
                                    'vmax_forward':[ODE,MCA,QSSA],
                                    'kcat_forward': [ODE, MCA, QSSA],
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


    def __init__(self, name, reactants, parameters=None, **kwargs):
        # FIXME dynamic linking, separaret parametrizations from model init
        # FIXME Reaction has a mechanism, and this is a mechanism
        KineticMechanism.__init__(self, name, reactants, parameters, **kwargs)

    def get_qssa_rate_expression(self):
        s = self.reactants.substrate.symbol
        p = self.reactants.product.symbol

        kms = self.parameters.km_substrate.symbol
        kmp = self.parameters.km_product.symbol

        keq = self.parameters.k_equilibrium.symbol
        if self.enzyme is not None:
            enzyme = self.reactants.enzyme.symbol
            kcat = self.parameters.kcat_forward.symbol
            vmaxf = kcat*enzyme
        else:
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
        rate_constants = [ 'k1_fwd','k1_bwd','k2_fwd','k2_bwd']
        if any([self.parameters[x].value is None for x in rate_constants] ):
            self.calculate_rate_constants()

        #TODO Better solution ?
        # Add enzyme and enzyme complex to the reactants of the mechanism
        self.reactants['enzyme'] = Reactant(self.name)
        self.reactants['enzyme_complex'] = Reactant('EC_'+self.name)

        s = self.reactants.substrate.symbol
        p = self.reactants.product.symbol
        e = self.reactants.enzyme.symbol
        es = self.reactants.enzyme_complex.symbol

        k1_fwd = self.parameters.k1_fwd.symbol
        k1_bwd = self.parameters.k1_bwd.symbol
        k2_fwd = self.parameters.k2_fwd.symbol
        k2_bwd = self.parameters.k2_bwd.symbol

        r1f = k1_fwd*e*s
        r1b = k1_bwd*es
        r2f = k2_fwd*es
        r2b = k2_bwd*e*p

        self.reaction_rates = TabDict([('r1f', r1f),
                                       ('r1b', r1b),
                                       ('r2f', r2f),
                                       ('r2b', r2b),
                                       ])

        self.expressions = {s: r1b - r1f,
                            p: r2f - r2b,
                            es: r1f - r1b - r2f + r2b,
                            e: r1b - r1f - r2b + r2f}

        parameters = [self.get_parameters_from_expression(expr)
                      for expr in self.expressions.values()]

        self.expression_parameters = set().union(*parameters)

    def calculate_rate_constants(self):

        params = self.parameters
        # Calcuate elementary rates
        #TODO Properly catch all possible parameter combinations that
        # parametrize the kinetics fully
        params.vmax_backward._symbol = params.vmax_forward.symbol \
                                       * params.km_product.symbol \
                                       / params.km_substrate.symbol \
                                       / params.k_equilibrium.symbol

        k1_bwd = params.vmax_backward.symbol / params.total_enzyme_concentration.symbol
        k2_fwd = params.vmax_forward.symbol  / params.total_enzyme_concentration.symbol

        k1_fwd = (k1_bwd + k2_fwd ) / params.km_substrate.symbol
        k2_bwd = (k1_bwd + k2_fwd ) / params.km_product.symbol

        self.parameters.k1_fwd._symbol = k1_fwd
        self.parameters.k1_bwd._symbol = k1_bwd
        self.parameters.k2_fwd._symbol = k2_fwd
        self.parameters.k2_bwd._symbol = k2_bwd

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
