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

    Rates = namedtuple('Rates',['k1_fwd',
                                'k1_bwd',
                                'k2_fwd',
                                'k2_bwd',
                                ])

    def __init__(self, name, substrates, parameters=None):
        # FIXME dynamic linking, separaret parametrizations from model init
        # FIXME Reaction has a mechanism, and this is a mechanism
        KineticMechanism.__init__(self,name, substrates, parameters)

    def calculate_rates(self):
        # Calcuate elementary rates

        # TODO: Check that the not None params are compatible
        # Param families/ set ?

        params = self.parameters
        k1_fwd = params.vmax_forward  / params.total_enzyme_concentration
        k2_fwd = params.vmax_backward / params.total_enzyme_concentration

        k1f = (k1_fwd + k2_fwd ) / params.km_substrate
        k2b = (k1_fwd + k2_fwd ) / params.km_product

        self.rates = self.Rates (k1_fwd = k1_fwd,
                                 k2_fwd = k2_fwd,
                                 k1_bwd = k1_bwd,
                                 k2_bwd = k2_bwd,
                                 )

        # self.set_dynamic_attribute_links(self._rates)


    def get_qssa_rate_expression(self):
        # FIXME dynamic linking, separated parametrizations from model init

        subs = self.substrates

        common_expression = sympify('1'                           \
                                    + '+'+subs.substrate+'/'      \
                                    +'km_substrate'+'_'+self.name \
                                    + '+'+subs.product+'/'        \
                                    +'km_product'+'_'+self.name  \
                                    )
        bwd_expression = sympify( 'vmax_backward'+'_'+self.name    \
                                  +'*'+subs.product               \
                                  +'/'+'km_product'+'_'+self.name)

        fwd_expression = sympify( 'vmax_forward'+'_'+self.name    \
                                  +'*'+subs.substrate             \
                                  +'/'+'km_substrate'+'_'+self.name)

        rate_expression = (fwd_expression-bwd_expression)/common_expression

        expressions = {subs.substrate: -rate_expression,
                       subs.product:    rate_expression}


        parameters = self.get_expression_parameters_from('parameters')

        return expressions, parameters


    def get_full_rate_expression(self):
        # Calculate rates uppon initialization
        if not hasattr(self,'rates'):
            self.calculate_rates()

        subs = self.substrates

        enzyme_complex = 'EC_'+self.name

        r1f = sympify(subs.substrates[0]+"*"+self.name+'*'+'k1_fwd'+self.name)
        r1b = sympify(enzyme_complex+'*k1_bwd'+self.name)
        r2f = sympify(enzyme_complex+'*k2_fwd'+self.name)
        r2b = sympify(subs.substrates[1]+"*"+self.name+'*k2_bwd'+self.name)

        expressions = {subs.substrate   : r1b - r1f,
                       subs.product     : r2f - r2b,
                       enzyme_complex   : r1f - r1b - r2f + r2b,
                       self.name        : r1b - r1f - r2b + r2f}


        parameters = self.get_expression_parameters_from('rates')

        return expressions,parameters




