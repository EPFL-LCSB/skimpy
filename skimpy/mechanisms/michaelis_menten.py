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
    Rates = namedtuple('Rates',['k1_fwd',
                                'k1_bwd',
                                'k2_fwd',
                                'k2_bwd',
                                ])

    def __init__(self, name, substrates, params):
        # FIXME dynamic linking, separaret parametrizations from model init
        # FIXME Reaction has a mechanism, and this is a mechanism
        KineticMechanism.__init__(name, substrates, params)
        # initialize from super class
        rates = {"k1f": k1f, "k1b": k1b, "k2f": k2f, "k2b": k2b}
        self.rates = rates

    def calculate_rates(self):
        k1b = self.params['v_max_ / self.params['E_tot']
        k2f = self.params['v_max_f'] / self.params['E_tot']
        k1f = (k1b + k2f) / self.params['K_S']
        k2b = (k1b + k2f) / self.params['K_P']

    def get_qssa_rate_expression(self):
        # FIXME dynamic linking, separated parametrizations from model init
        denominator = sympify('1+' + self.substrates[1] + '/K_P_' + self.name \
                              + "+" + self.substrates[0] + '/K_S_' + self.name)
        reverse_flux = sympify(
            'v_max_r_' + self.name + "*" + self.substrates[1] \
            + "/K_P_" + self.name)
        forward_flux = sympify(
            'v_max_f_' + self.name + "*" + self.substrates[0] \
            + "/K_S_" + self.name)

        rate = (forward_flux - reverse_flux) / denominator

        expressions = {self.substrates[0]: (-1.0) * rate,
                       self.substrates[1]: rate}

        variables = [self.substrates[0], self.substrates[1]]

        parameters = {}
        for this_key in self.params:
            parameters[this_key + '_' + self.name] = self.params[this_key]

        return variables, expressions, parameters


    def get_full_rate_expression(self):

        enzyme_complex = 'EC_'+self.name

        r1f = sympify(self.substrates[0]+"*"+self.name+'*'+'k1f_'+self.name)
        r1b = sympify(enzyme_complex+'*k1b_'+self.name)
        r2f = sympify(enzyme_complex+'*k2f_'+self.name)
        r2b = sympify(self.substrates[1]+"*"+self.name+'*k2b_'+self.name)

        expressions = {self.substrates[0] : r1b - r1f,
                      self.substrates[1]  : r2f - r2b,
                      enzyme_complex : r1f - r1b - r2f + r2b,
                      self.name      : r1b - r1f - r2b + r2f}
        variables   = [self.substrates[0], self.substrates[1] , enzyme_complex, self.name]

        parameters = {}
        for this_key in self.rates:
            parameters[this_key+'_'+self.name] = self.rates[this_key]

        return variables,expressions,parameters




