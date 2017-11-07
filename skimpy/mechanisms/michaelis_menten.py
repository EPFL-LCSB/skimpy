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
from .mechanism import KineticMechanism
from ..core.reactions import Reaction
from ..utils.tabdict import TabDict
from collections import namedtuple


class ReversibleMichaelisMenten(KineticMechanism):
    """A reversible uni uni reaction enzmye class"""

    Substrates = namedtuple('Substrates', ['substrate', 'product'])

    Parameters = namedtuple('Parameters', ['vmax_forward',
                                           'vmax_backward',
                                           'km_substrate',
                                           'km_product',
                                           'thermo_displacement',
                                           'k_equilibrium',
                                           'total_enzyme_concentration',
                                           ])
    Parameters.__new__.__defaults__ = (None,) * len(Parameters._fields)

    RateConstants = namedtuple('RateConstants',['k1_fwd',
                                                'k1_bwd',
                                                'k2_fwd',
                                                'k2_bwd',
                                                ])

    def __init__(self, name, substrates, parameters=None):
        # FIXME dynamic linking, separaret parametrizations from model init
        # FIXME Reaction has a mechanism, and this is a mechanism
        KineticMechanism.__init__(self,name, substrates, parameters)

    def get_qssa_rate_expression(self):
        subs = self.substrates

        common_denominator = sympify('1'                           \
                                    + '+'+subs.substrate+'/'      \
                                    +'km_substrate'+'_'+self.name \
                                    + '+'+subs.product+'/'        \
                                    +'km_product'+'_'+self.name  \
                                    )
        bwd_nominator = sympify( 'vmax_backward'+'_'+self.name    \
                                  +'*'+subs.product               \
                                  +'/'+'km_product'+'_'+self.name)

        fwd_nominator = sympify( 'vmax_forward'+'_'+self.name    \
                                  +'*'+subs.substrate             \
                                  +'/'+'km_substrate'+'_'+self.name)

        forward_rate_expression = fwd_nominator/common_denominator
        backward_rate_expression = bwd_nominator/common_denominator
        rate_expression = forward_rate_expression-backward_rate_expression

        self.reaction_rates = TabDict([('v_net',rate_expression),
                                       ('v_fwd', forward_rate_expression),
                                       ('v_bwd', backward_rate_expression),
                                       ])

        expressions = {subs.substrate: -rate_expression,
                       subs.product:    rate_expression}


        parameters = self.get_expression_parameters_from('parameters')

        return expressions, parameters


    def get_full_rate_expression(self):
        # Calculate rates uppon initialization
        if not hasattr(self,'rate_constants'):
            self.calculate_rate_constants()

        subs = self.substrates

        enzyme_complex = 'EC_'+self.name

        r1f = sympify(subs.substrate+"*"+self.name+'*'+'k1_fwd_'+self.name)
        r1b = sympify(enzyme_complex+'*k1_bwd_'+self.name)
        r2f = sympify(enzyme_complex+'*k2_fwd_'+self.name)
        r2b = sympify(subs.product  +"*"+self.name+'*k2_bwd_'+self.name)

        self.reaction_rates = TabDict([('r1f',r1f),
                                       ('r1b', r1b),
                                       ('r2f', r2f),
                                       ('r2b', r2b),
                                       ])

        expressions = {subs.substrate: r1b - r1f,
                       subs.product: r2f - r2b,
                       enzyme_complex: r1f - r1b - r2f + r2b,
                       self.name: r1b - r1f - r2b + r2f}

        parameters = self.get_expression_parameters_from('rate_constants')
        return expressions, parameters


    def calculate_rate_constants(self):
        # Calcuate elementary rates

        # TODO: Check that the not None params are compatible
        # Param families/ set ?

        params = self.parameters
        k1_bwd = params.vmax_backward / params.total_enzyme_concentration
        k2_fwd = params.vmax_forward  / params.total_enzyme_concentration

        k1_fwd = (k1_bwd + k2_fwd ) / params.km_substrate
        k2_bwd = (k1_bwd + k2_fwd ) / params.km_product

        self.rate_constants = self.RateConstants (k1_fwd = k1_fwd,
                                                  k2_fwd = k2_fwd,
                                                  k1_bwd = k1_bwd,
                                                  k2_bwd = k2_bwd,
                                                  )

        # self.set_dynamic_attribute_links(self._rates)
