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


class RandBiBiReversibleMichaelisMenten(KineticMechanism):
    """A reversible random ordered bi bi reaction enzmye class"""

    Substrates = namedtuple('Substrates', ['substrate1',
                                           'substrate2',
                                           'product1',
                                           'product2'])

    Parameters = namedtuple('Parameters', ['vmax_forward',
                                           'k_equilibrium',
                                           'kmi_substrate1',
                                           'kmi_substrate2',
                                           'km_substrate2',
                                           'kmi_product1',
                                           'kmi_product2',
                                           'km_product2',
                                           #'vmax_backward',
                                           #'total_enzyme_concentration',
                                           ])

    Parameters.__new__.__defaults__ = (None,) * len(Parameters._fields)

    RateConstants = namedtuple('RateConstants',[])
    ElementaryReactions = namedtuple('ElementaryReactions',[])


    def __init__(self, name, substrates, parameters=None):
        # FIXME dynamic linking, separaret parametrizations from model init
        # FIXME Reaction has a mechanism, and this is a mechanism
        KineticMechanism.__init__(self,name, substrates, parameters)

    def get_qssa_rate_expression(self):
        subs = self.substrates

        common_denominator = sympify('1'
                                      + '+'+subs.substrate1+'/'
                                      + 'kmi_substrate1'+'_'+self.name
                                      + '+' + subs.substrate2 + '/'
                                      + 'kmi_substrate2' + '_' + self.name
                                      + '+' + subs.product1 + '/'
                                      + 'kmi_product1' + '_' + self.name
                                      + '+'+subs.product2+'/'
                                      + 'kmi_product2'+'_'+self.name
                                      + '+' + subs.substrate1
                                      + '*' + subs.substrate2
                                      + '/' + 'kmi_substrate1' + '_' + self.name
                                      + '/' + 'km_substrate2' + '_' + self.name
                                      + '+' + subs.product1
                                      + '*' + subs.product2
                                      + '/' + 'kmi_product1' + '_' + self.name
                                      + '/' + 'km_product2' + '_' + self.name
                                     )

        bwd_nominator = sympify( 'vmax_forward'+'_'+self.name
                                  + '/k_equilibrium'+'_'+self.name
                                  + '*' + subs.product1
                                  + '*' + subs.product2
                                  + '/' + 'kmi_substrate1' + '_' + self.name
                                  + '/' + 'km_substrate2' + '_' + self.name
                                 )

        fwd_nominator = sympify( 'vmax_forward'+'_'+self.name
                                  + '*' + subs.substrate1
                                  + '*' + subs.substrate2
                                  + '/' + 'kmi_substrate1' + '_'+self.name
                                  + '/' + 'km_substrate2' + '_' + self.name
                                )

        forward_rate_expression = fwd_nominator/common_denominator
        backward_rate_expression = bwd_nominator/common_denominator
        rate_expression = forward_rate_expression-backward_rate_expression

        self.reaction_rates = TabDict([('v_net', rate_expression),
                                       ('v_fwd', forward_rate_expression),
                                       ('v_bwd', backward_rate_expression),
                                       ])

        expressions = {subs.substrate1: -rate_expression,
                       subs.substrate2: -rate_expression,
                       subs.product1:    rate_expression,
                       subs.product2:    rate_expression
                       }

        parameters = self.get_expression_parameters_from('parameters')

        return expressions, parameters


    def get_full_rate_expression(self):
        raise NotImplementedError

    def calculate_rate_constants(self):
        raise NotImplementedError
