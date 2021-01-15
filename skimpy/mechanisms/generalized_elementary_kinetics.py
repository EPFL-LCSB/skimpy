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

from sympy import Symbol, exp
from .mechanism import KineticMechanism
from ..utils.tabdict import TabDict
from collections import namedtuple
from ..core.itemsets import make_parameter_set, make_reactant_set
from skimpy.utils.general import make_subclasses_dict
from ..utils.namespace import *
from .utils import stringify_stoichiometry


def make_generalized_elementary_kinetics(stoichiometry, generalized_reactants):
    """
    Creates a reversible N-M GeneralizedElementaryKinetics class

    GeneralizedElementaryKinetics iplementing the kinetic scheme to account for non ideal effects in large scale kinetic model:

    Weilandt, Daniel R., and Vassily Hatzimanikatis. "Particle-based simulation reveals macromolecular
    crowding effects on the Michaelis-Menten mechanism." Biophysical Journal (2019).
    
    :param stoichiometry is a list of the reaction stoichioemtry
    :param generalized_reactants is a list of reactants that contribute in a powerlaw form:
           v = v_0 * exp(beta_forward)
                ... * (x_1/x_1_0)^(alpha_f_x_1 )
                ... * (x_n/x_n_0)^(alpha_f_x_n )
        with v_0 beeing the reversible Massaction flux 

    """
    # This refresh all subclasses and fetches already create mechanism classes
    ALL_MECHANISM_SUBCLASSES = make_subclasses_dict(KineticMechanism)

    new_class_name = "GEEK" \
                     + "_{0}".format(stringify_stoichiometry(stoichiometry))

    if new_class_name in ALL_MECHANISM_SUBCLASSES.keys():
        return ALL_MECHANISM_SUBCLASSES[new_class_name]

    class GeneralizedElementaryKinetics(KineticMechanism):
        """
        A reversible N-M GEEK class
        Implementing the kinetic scheme to account for non ideal effects in large scale kinetic model
        Weilandt, Daniel R., and Vassily Hatzimanikatis. "Particle-based simulation reveals macromolecular
        crowding effects on the Michaelis-Menten mechanism." Biophysical Journal (2019).

        """

        suffix = "_{0}_{1}".format(stringify_stoichiometry(stoichiometry), generalized_reactants)

        reactant_list = []
        parameter_list = {'vmax_forward': [ODE, MCA, QSSA],
                          'k_equilibrium': [ODE, MCA, QSSA],
                          'beta_forward': [ODE, MCA, QSSA],
                          'beta_reverse': [ODE, MCA, QSSA],}

        reactant_parameter_links = {}
        reactant_stoichiometry = {}

        num_substrates = 1
        num_products = 1
        for s in stoichiometry:
            if s < 0:
                substrate = 'substrate{}'.format(num_substrates)
                reactant_list.append(substrate)
                reactant_stoichiometry[substrate] = float(s)
                num_substrates += 1

            if s > 0:
                product = 'product{}'.format(num_products)
                reactant_list.append(product)
                reactant_stoichiometry[product] = float(s)
                num_products += 1
           

        for gr in generalized_reactants:

            alpha_f = 'alpha_forward_{}'.format(gr)
            parameter_list[alpha_f] =  [ODE, MCA, QSSA]

            alpha_r = 'alpha_reverse_{}'.format(gr)
            parameter_list[alpha_r] =  [ODE, MCA, QSSA]

            X_0 = '{}_0'.format(gr)
            parameter_list[X_0] = [ODE, MCA, QSSA]  # TODO Does this work if multiple times defined

            reactant_parameter_links[gr] = (alpha_f, alpha_r, X_0)


        Reactants = make_reactant_set(__name__ + suffix, reactant_list)

        Parameters = make_parameter_set(__name__ + suffix, parameter_list)

        Inhibitors = make_reactant_set(__name__ + suffix, generalized_reactants)

        #ElementaryReactions = namedtuple('ElementaryReactions', [])

        def __init__(self, name, reactants, parameters=None, inhibitors=None, **kwargs):
            # FIXME dynamic linking, separaret parametrizations from model init
            # FIXME Reaction has a mechanism, and this is a mechanism
            KineticMechanism.__init__(self, name, reactants, parameters, **kwargs)

            if inhibitors is not None:
                self.inhibitors = inhibitors
            else:
                self.inhibitors = self.Inhibitors(**dict(zip(generalized_reactants, generalized_reactants)))

        def get_qssa_rate_expression(self):

            substrates = {k: r for k, r in self.reactants.items()
                          if k.startswith('substrate')}

            products = {k: r for k, r in self.reactants.items()
                        if k.startswith('product')}

            kf = self.parameters.vmax_forward.symbol
            Keq = self.parameters.k_equilibrium.symbol

            beta_f = self.parameters.beta_forward.symbol
            beta_r = self.parameters.beta_reverse.symbol


            forward_rate_expression = kf * exp(beta_f)
            backward_rate_expression = kf / Keq * exp(beta_r)

            for this_type, this_substrate in substrates.items():
                s = this_substrate.symbol
                forward_rate_expression *= s


            for this_type, this_product in products.items():
                p = this_product.symbol
                backward_rate_expression *= p

            for this_gr in self.inhibitors.values():

                alpha_f, alpha_r, gr_0 = [self.parameters[p].symbol for p in
                                          self.reactant_parameter_links[this_gr.name]]

                # FIXME  get this from a proper defined set
                gr = this_gr.symbol
                forward_rate_expression *= (gr/gr_0)**alpha_f
                backward_rate_expression *= (gr/gr_0)**alpha_r


            rate_expression = forward_rate_expression - backward_rate_expression

            self.reaction_rates = TabDict([('v_net', rate_expression),
                                           ('v_fwd', forward_rate_expression),
                                           ('v_bwd', backward_rate_expression),
                                           ])

            expressions = {}

            for type, this_substrate in substrates.items():
                s = this_substrate.symbol
                stoich = self.reactant_stoichiometry[type]
                expressions[s] = stoich * rate_expression

            for type, this_product in products.items():
                p = this_product.symbol
                stoich = self.reactant_stoichiometry[type]
                expressions[p] = stoich * rate_expression

            self.expressions = expressions
            self.expression_parameters = self.get_parameters_from_expression(rate_expression)

        def update_qssa_rate_expression(self):

            substrates = {k: r for k, r in self.reactants.items()
                          if k.startswith('substrate')}

            products = {k: r for k, r in self.reactants.items()
                        if k.startswith('product')}

            for type, this_substrate in substrates.items():
                s = this_substrate.symbol
                stoich = self.reactant_stoichiometry[type]
                self.expressions[s] = stoich * self.reaction_rates['v_net']

            for type, this_product in products.items():
                p = this_product.symbol
                stoich = self.reactant_stoichiometry[type]
                self.expressions[p] = stoich * self.reaction_rates['v_net']

        """"
        Not implemented not necessary for elementary mechanisms
        """

        def get_full_rate_expression(self):
            raise NotImplementedError

        def calculate_rate_constants(self):
            raise NotImplementedError

    GeneralizedElementaryKinetics.__name__ += GeneralizedElementaryKinetics.suffix

    return GeneralizedElementaryKinetics

