"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2021 Laboratory of Computational Systems Biotechnology (LCSB),
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


class SimpleRegulatedGeneExpression(KineticMechanism):
    """
    Apprxoimate gene regulation according to XX et al
    e.g.: A + B -> C

    """

    reactant_list = ['product',]
    Reactants = make_reactant_set(__name__, reactant_list)

    # Todo rename to regulators!
    regulator_list = ['factor_1', ]
    SingleRegulator = make_reactant_set(__name__, regulator_list)

    regulator_list = ['factor_1', 'factor_2']
    Regulators = make_reactant_set(__name__, regulator_list)

    Parameters = make_parameter_set(__name__,
                                    {
                                        'basal_transcription_rate':[ODE,MCA,QSSA],
                                        'transcription_rate_factor_1':[ODE,MCA,QSSA],
                                        'transcription_rate_factor_2':[ODE,MCA,QSSA],
                                        'transcription_rate_factor_1_factor_2':[ODE,MCA,QSSA],
                                        'hill_coefficient_factor_1': [ODE, MCA, QSSA],
                                        'hill_coefficient_factor_2': [ODE, MCA, QSSA],
                                        'kd_factor_1': [ODE, MCA, QSSA],
                                        'kd_factor_2': [ODE, MCA, QSSA],
                                        'cooperativity':[ODE,MCA,QSSA],
                                    })

    reactant_stoichiometry = {'product': 1
                              }

    parameter_reactant_links = {}


    ElementaryReactions = namedtuple('ElementaryReactions',[])

    def __init__(self, name, reactants, inhibitors, parameters=None, **kwargs):
        KineticMechanism.__init__(self, name, reactants, parameters, inhibitors=inhibitors, **kwargs)

    def get_qssa_rate_expression(self):
        products = TabDict([(k, self.reactants[k])
                            for k in self.reactant_list
                            if k.startswith('product')])

        basal_transcription_rate = self.parameters.basal_transcription_rate.symbol

        if not hasattr(self.inhibitors,'factor_1'):
            factor_1 = 0
            transcription_rate_factor_1 = 0
            hill_coefficient_factor_1 = 1
            kd_factor_1 = 1
        else:
            factor_1 = self.inhibitors.factor_1.symbol
            transcription_rate_factor_1 = \
                self.parameters.transcription_rate_factor_1.symbol
            hill_coefficient_factor_1 = \
                self.parameters.hill_coefficient_factor_1.symbol
            kd_factor_1 = \
                self.parameters.kd_factor_1.symbol

        if not hasattr(self.inhibitors,'factor_2'):
            factor_2 = 0
            transcription_rate_factor_2 = 0
            hill_coefficient_factor_2 = 1
            transcription_rate_factor_1_factor_2 = 0
            kd_factor_2 = 1
            cooperativity = 0
        else:
            factor_2 = self.inhibitors.factor_2.symbol
            transcription_rate_factor_2 = \
                self.parameters.transcription_rate_factor_2.symbol
            hill_coefficient_factor_2 = \
                self.parameters.hill_coefficient_factor_2.symbol
            transcription_rate_factor_1_factor_2 = \
                self.parameters.transcription_rate_factor_1_factor_2.symbol
            cooperativity = self.parameters.cooperativity.symbol
            kd_factor_2 = \
                self.parameters.kd_factor_2.symbol

        nominator = basal_transcription_rate + \
                    transcription_rate_factor_1*(factor_1/kd_factor_1)**hill_coefficient_factor_1 + \
                    transcription_rate_factor_2*(factor_2/kd_factor_2)**hill_coefficient_factor_2 + \
                    transcription_rate_factor_1_factor_2*cooperativity *\
                    (factor_1/kd_factor_1)**hill_coefficient_factor_1 * \
                    (factor_2/kd_factor_2)**hill_coefficient_factor_2

        denominator = 1 + (factor_1/kd_factor_1)**hill_coefficient_factor_1 \
                      + (factor_1/kd_factor_1)**hill_coefficient_factor_1 + \
                      cooperativity * (factor_1/kd_factor_1)**hill_coefficient_factor_1 * \
                      (factor_2/kd_factor_2)**hill_coefficient_factor_2


        forward_rate_expression = nominator/denominator

        backward_rate_expression = 0.0
        rate_expression = forward_rate_expression-backward_rate_expression

        self.reaction_rates = TabDict([('v_net', rate_expression),
                                       ('v_fwd', forward_rate_expression),
                                       ('v_bwd', backward_rate_expression),
                                       ])

        expressions = {}

        # for type, this_substrate in substrates.items():
        #     s = this_substrate.symbol
        #     stoich = self.reactant_stoichiometry[type]
        #     expressions[s] = stoich*rate_expression

        for type, this_product in products.items():
            p = this_product.symbol
            stoich = self.reactant_stoichiometry[type]
            expressions[p] = stoich * rate_expression

        self.expressions = expressions
        self.expression_parameters = self.get_parameters_from_expression(rate_expression)

    def update_qssa_rate_expression(self):

        # substrates = {k:r for k,r in self.reactants.items()
        #               if k.startswith('substrate')}

        products= {k:r for k,r in self.reactants.items()
                   if k.startswith('product')}

        # for type, this_substrate in substrates.items():
        #     s = this_substrate.symbol
        #     stoich = self.reactant_stoichiometry[type]
        #     self.expressions[s] = stoich*self.reaction_rates['v_net']

        for type, this_product in products.items():
            p = this_product.symbol
            stoich = self.reactant_stoichiometry[type]
            self.expressions[p] = stoich*self.reaction_rates['v_net']


    """"
    This kinetics has no detailed mechanism 
    """
    def get_full_rate_expression(self):
        raise NotImplementedError

    def calculate_rate_constants(self):
        raise NotImplementedError

