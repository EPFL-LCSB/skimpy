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
from ..utils.general import check_is_symbol
from ..mechanisms.mechanism import KineticMechanism
from ..core.itemsets import make_parameter_set, make_reactant_set
from ..utils.namespace import *


class ExpressionModifier(object):
    """
    This class describes a modifier to an expression, like a boundary condition
    or constraint.
    For example, changing a rate to a constant (boundary condition), or linking
    it to another variable of the model (constraint).
    It accepts as an argument a modifier.

    A modifier is a function which will look at all your expressions, and
    apply its transformation to them. As a result, its arguments have to be a
    TabDict of expressions, such as KinModel.ODEFun.expressions
    """

    prefix = 'MOD'

    def __init__(self, name, reaction=None, modifier=None):
        self._name = name
        if modifier is not None:
            self._modifier = modifier

        self._reaction = reaction

    def __call__(self,expressions):
        self.modifier(expressions)

    @property
    def modifier(self):
        return self._modifier

    def link(self,model):
        """
        Link the modifier to a model, to gain awareness of the inner/outer
        variables
        :param model:
        :return:
        """

        self.model = model

    @property
    def name(self):
        return self.prefix +'_' + self._name

    @name.setter
    def name(self, value):
        if value.startswith(self.prefix):
            value = value[len(self.prefix):]
        self._name = value


class BoundaryCondition(ExpressionModifier):
    """
    We differentiate boundary conditions as modifiers that define the boundaries
    of the observed system.
    """

    prefix = 'BC'

    def __init__(self, name, modifier = None, reaction=None):
        ExpressionModifier.__init__(self, name, modifier)


class ConstantConcentration(BoundaryCondition):
    """

    """

    prefix = 'CC'

    def __init__(self, reactant, name = None, reaction=None):

        # Is the reactant constant it is not a variable anymore
        if name is None:
            name = reactant.name

        BoundaryCondition.__init__(self, name = name)

        # Modify the reactant
        reactant.type = PARAMETER
        self.reactant = reactant

    def modifier(self, expressions):
        """
        Set the rate to 0
        :param expressions:
        :return:
        """
        expressions[self.reactant.symbol] = expressions[self.reactant.symbol] * 0.0

    def __del__(self):
        self.reactant.type = VARIABLE


class AdditiveConcentrationRate(ExpressionModifier):
    """
    Add a concentration rate term to your rate expression
    """
    # FIXME Please give us an alternate name we _REALLY_ don't like it

    prefix = 'ADDCR'

    def __init__(self, reactant, flux_value, name=None, reaction=None):

        if name is None:
            name = reactant.__str__()

        ExpressionModifier.__init__(self, name=name)

        self.reactant = reactant
        self.flux_value = flux_value

    def modifier(self, expressions):
        """
        Add to the rate expression
        :param expressions:
        :return:
        """
        sym_value = sympify(self.flux_value)
        expressions[self.reactant.symbol] = expressions[self.reactant.symbol] + sym_value


class BoundaryFlux(BoundaryCondition,AdditiveConcentrationRate):

    prefix = "BF"

    def __init__(self, reactant, flux_value, reaction=None):
        # TODO: Find a way to make sure the flux_value does not depend on an
        # inner variable
        self.check_dependency(flux_value)
        AdditiveConcentrationRate.__init__(self, reactant, flux_value)

    def check_dependency(self, expression):
        # TODO: Implement
        pass


"""
Reaction modifiers 
"""

class FirstOrderSmallMoleculeModifier(KineticMechanism,ExpressionModifier):

    prefix = "HSM"

    Reactants = make_reactant_set(__name__, ['small_molecule'])

    Parameters = make_parameter_set(    __name__,
                                        { })

    parameter_reactant_links = {}

    def __init__(self, small_molecule, mechanism_stoichiometry, name=None, reaction=None):

        if name is None:
            name = small_molecule.__repr__()

        reactants = self.Reactants(small_molecule=small_molecule)
        parameters = self.Parameters()
        KineticMechanism.__init__(self, name, reactants, parameters)

        if type(mechanism_stoichiometry) is dict:
            self.reactant_stoichiometry = mechanism_stoichiometry
        else:
            self.reactant_stoichiometry = {'small_molecule':
                                               float(mechanism_stoichiometry)}

    def modifier(self, expressions):
        """
        change the flux reaction rate expressions
        :param expression: {vnet, vfwd, vbwd}
        :return:
        """
        # First oder modification of the of Keq
        # expressions = TabDict([('v_net', rate_expression),
        #                        ('v_fwd', forward_rate_expression),
        #                        ('v_bwd', backward_rate_expression),
        #                        ])

        if self.reactant_stoichiometry['small_molecule'] < 0:
            expressions['v_fwd'] = expressions['v_fwd']\
                                   * self.get_qssa_rate_expression()**-self.reactant_stoichiometry['small_molecule']

        if self.reactant_stoichiometry['small_molecule'] > 0:
            expressions['v_bwd'] = expressions['v_bwd'] \
                                   * self.get_qssa_rate_expression()**self.reactant_stoichiometry['small_molecule']

        expressions['v_net'] = expressions['v_fwd'] - expressions['v_bwd']

    def get_qssa_rate_expression(self):
        sm = self.reactants.small_molecule.symbol
        return sm

    def update_qssa_rate_expression(self):
        return None


    def get_full_rate_expression(self):
        raise NotImplementedError

    def calculate_rate_constants(self):
        raise NotImplementedError


class DisplacementSmallMoleculeModifier(KineticMechanism,ExpressionModifier):

    prefix = "DSM"

    Reactants = make_reactant_set(__name__, ['small_molecule',])

    Parameters = make_parameter_set(    __name__,
                                        { })

    parameter_reactant_links = {}

    def __init__(self, small_molecule, mechanism_stoichiometry, name=None, reaction=None):

        if name is None:
            name = small_molecule.__str__()

        reactants = self.Reactants(small_molecule=small_molecule,)
        parameters = self.Parameters()
        KineticMechanism.__init__(self, name, reactants, parameters)

        # TODO Unify between skimpy versions
        if type(mechanism_stoichiometry) is dict:
            self.reactant_stoichiometry = mechanism_stoichiometry
        else:
            self.reactant_stoichiometry = {'small_molecule':
                                               float(mechanism_stoichiometry)}

    def modifier(self, expressions):
        """
        change the flux reaction rate expressions
        :param expression: {vnet, vfwd, vbwd}
        :return:
        """
        # Modification of the of Keq
        # expressions = TabDict([('v_net', rate_expression),
        #                        ('v_fwd', forward_rate_expression),
        #                        ('v_bwd', backward_rate_expression),
        #                        ])

        expressions['v_bwd'] = expressions['v_bwd'] \
                                * self.get_qssa_rate_expression()**self.reactant_stoichiometry['small_molecule']

        expressions['v_net'] = expressions['v_fwd'] - expressions['v_bwd']

    def get_qssa_rate_expression(self):
        sm = self.reactants.small_molecule.symbol
        return sm

    def update_qssa_rate_expression(self):
        return None

    def get_full_rate_expression(self):
        raise NotImplementedError

    def calculate_rate_constants(self):
        raise NotImplementedError

"""
Activators and inhibitors
"""

class ActivationModifier(KineticMechanism,ExpressionModifier):

    prefix = "AM"

    Reactants = make_reactant_set(__name__, ['activator',])

    Parameters = make_parameter_set(__name__, {'k_activation': [ODE, MCA, QSSA],})

    parameter_reactant_links = {'k_activation':'activator'}

    def __init__(self, activator, name=None, k_activation=None, reaction=None):

        if name is None:
            name = self.prefix+'_'+activator.__str__()

        if reaction is None:
            suffix = name
        else:
            suffix = name+'_'+reaction.name

        reactants = self.Reactants(activator=activator,)
        parameters = self.Parameters(k_activation=k_activation)

        for name, p in parameters.items():
            p.suffix = suffix

        KineticMechanism.__init__(self, name, reactants, parameters)

        self.link_parameters_and_reactants()

        self.reactant_stoichiometry = {'activator': 0 }

    def modifier(self, expressions):
        """
        change the flux reaction rate expressions
        :param expression: {vnet, vfwd, vbwd}
        :return:
        """
        # Modification of the of Keq
        # expressions = TabDict([('v_net', rate_expression),
        #                        ('v_fwd', forward_rate_expression),
        #                        ('v_bwd', backward_rate_expression),
        #                        ])
        activation = 1 + self.get_qssa_rate_expression()
        expressions['v_bwd'] = expressions['v_bwd'] * activation
        expressions['v_fwd'] = expressions['v_fwd'] * activation
        expressions['v_net'] = expressions['v_fwd'] - expressions['v_bwd']

    def get_qssa_rate_expression(self):
        a = self.reactants.activator.symbol
        k = self.parameters.k_activation.symbol
        return a/k

    def update_qssa_rate_expression(self):
        return None

    def get_full_rate_expression(self):
        raise NotImplementedError

    def calculate_rate_constants(self):
        raise NotImplementedError


class InhibitionModifier(KineticMechanism,ExpressionModifier):

    prefix = "IM"

    Reactants = make_reactant_set(__name__, ['inhibitor',])

    Parameters = make_parameter_set(__name__, {'k_inhibition': [ODE, MCA, QSSA],})

    parameter_reactant_links = {'k_inhibition':'inhibitor'}

    def __init__(self, inhibitor, name=None, k_inhibition=None, reaction=None):

        if name is None:
            name = self.prefix+'_'+inhibitor.__str__()

        if reaction is None:
            suffix = name
        else:
            suffix = name+'_'+reaction.name

        reactants = self.Reactants(inhibitor=inhibitor,)
        parameters = self.Parameters(k_inhibition=k_inhibition)

        for name, p in parameters.items():
            p.suffix = suffix

        KineticMechanism.__init__(self, name, reactants, parameters)
        self.link_parameters_and_reactants()

        self.reactant_stoichiometry = {'inhibitor': 0 }

    def modifier(self, expressions):
        """
        change the flux reaction rate expressions
        :param expression: {vnet, vfwd, vbwd}
        :return:
        """
        # Modification of the of Keq
        # expressions = TabDict([('v_net', rate_expression),
        #                        ('v_fwd', forward_rate_expression),
        #                        ('v_bwd', backward_rate_expression),
        #                        ])
        inhibition = 1 + self.get_qssa_rate_expression()
        expressions['v_bwd'] = expressions['v_bwd'] / inhibition
        expressions['v_fwd'] = expressions['v_fwd'] / inhibition
        expressions['v_net'] = expressions['v_fwd'] - expressions['v_bwd']

    def get_qssa_rate_expression(self):
        a = self.reactants.inhibitor.symbol
        k = self.parameters.k_inhibition.symbol
        return a/k

    def update_qssa_rate_expression(self):
        return None

    def get_full_rate_expression(self):
        raise NotImplementedError

    def calculate_rate_constants(self):
        raise NotImplementedError


class HillActivationModifier(KineticMechanism,ExpressionModifier):

    prefix = "HAM"

    Reactants = make_reactant_set(__name__, ['activator',])

    Parameters = make_parameter_set(__name__, {'k_activation': [ODE, MCA, QSSA],
                                               'a_max': [ODE, MCA, QSSA],
                                                'hill_coefficient':[ODE, MCA, QSSA], })

    parameter_reactant_links = {'k_activation':'activator'}

    def __init__(self, activator, name=None,
                 k_activation=None,
                 a_max=None,
                 hill_coefficient=None,
                 reaction=None):

        if name is None:
            name = self.prefix+'_'+activator.__str__()

        if reaction is None:
            suffix = name
        else:
            suffix = name+'_'+reaction.name

        reactants = self.Reactants(activator=activator,)
        parameters = self.Parameters(k_activation=k_activation,
                                     a_max=a_max,
                                     hill_coefficient=hill_coefficient)

        for name, p in parameters.items():
            p.suffix = suffix

        KineticMechanism.__init__(self, name, reactants, parameters)

        self.link_parameters_and_reactants()

        self.reactant_stoichiometry = {'activator': 0 }

    def modifier(self, expressions):
        """
        change the flux reaction rate expressions
        :param expression: {vnet, vfwd, vbwd}
        :return:
        """
        # Modification of the of Keq
        # expressions = TabDict([('v_net', rate_expression),
        #                        ('v_fwd', forward_rate_expression),
        #                        ('v_bwd', backward_rate_expression),
        #                        ])
        activation = 1 + self.get_qssa_rate_expression()
        expressions['v_bwd'] = expressions['v_bwd'] * activation
        expressions['v_fwd'] = expressions['v_fwd'] * activation
        expressions['v_net'] = expressions['v_fwd'] - expressions['v_bwd']

    def get_qssa_rate_expression(self):
        a = self.reactants.activator.symbol
        k = self.parameters.k_activation.symbol
        a_max = self.parameters.a_max.symbol
        h = self.parameters.hill_coefficient.symbol
        return a_max * ( (a/k)**h / (1 + (a/k)**h ))

    def update_qssa_rate_expression(self):
        return None

    def get_full_rate_expression(self):
        raise NotImplementedError

    def calculate_rate_constants(self):
        raise NotImplementedError


class SimpleHillActivationModifier(KineticMechanism,ExpressionModifier):

    prefix = "SHAM"

    Reactants = make_reactant_set(__name__, ['activator',])

    Parameters = make_parameter_set(__name__, {'k_activation': [ODE, MCA, QSSA],
                                                'hill_coefficient':[ODE, MCA, QSSA], })

    parameter_reactant_links = {'k_activation':'activator'}

    def __init__(self, activator, name=None,
                 k_activation=None,
                 a_max=None,
                 hill_coefficient=None,
                 reaction=None):

        if name is None:
            name = self.prefix+'_'+activator.__str__()

        if reaction is None:
            suffix = name
        else:
            suffix = name+'_'+reaction.name

        reactants = self.Reactants(activator=activator,)
        parameters = self.Parameters(k_activation=k_activation,
                                     hill_coefficient=hill_coefficient)

        for name, p in parameters.items():
            p.suffix = suffix

        KineticMechanism.__init__(self, name, reactants, parameters)

        self.link_parameters_and_reactants()

        self.reactant_stoichiometry = {'activator': 0 }

    def modifier(self, expressions):
        """
        change the flux reaction rate expressions
        :param expression: {vnet, vfwd, vbwd}
        :return:
        """
        # Modification of the of Keq
        # expressions = TabDict([('v_net', rate_expression),
        #                        ('v_fwd', forward_rate_expression),
        #                        ('v_bwd', backward_rate_expression),
        #                        ])
        activation = 1 + self.get_qssa_rate_expression()
        expressions['v_bwd'] = expressions['v_bwd'] * activation
        expressions['v_fwd'] = expressions['v_fwd'] * activation
        expressions['v_net'] = expressions['v_fwd'] - expressions['v_bwd']

    def get_qssa_rate_expression(self):
        a = self.reactants.activator.symbol
        k = self.parameters.k_activation.symbol
        h = self.parameters.hill_coefficient.symbol
        return (a/k)**h

    def update_qssa_rate_expression(self):
        return None

    def get_full_rate_expression(self):
        raise NotImplementedError

    def calculate_rate_constants(self):
        raise NotImplementedError



class HillInhibitionModifier(KineticMechanism,ExpressionModifier):

    prefix = "HIM"

    Reactants = make_reactant_set(__name__, ['inhibitor',])

    Parameters = make_parameter_set(__name__, {'k_inhibition': [ODE, MCA, QSSA],
                                               'hill_coefficient':[ODE, MCA, QSSA], })

    parameter_reactant_links = {'k_inhibition':'inhibitor'}

    def __init__(self, inhibitor, name=None,
                 k_inhibition=None,
                 hill_coefficient=None,
                 reaction=None):

        if name is None:
            name = self.prefix+'_'+inhibitor.__str__()

        if reaction is None:
            suffix = name
        else:
            suffix = name+'_'+reaction.name

        reactants = self.Reactants(inhibitor=inhibitor,)
        parameters = self.Parameters(k_inhibition=k_inhibition,
                                     hill_coefficient=hill_coefficient)

        for name, p in parameters.items():
            p.suffix = suffix

        KineticMechanism.__init__(self, name, reactants, parameters)

        self.link_parameters_and_reactants()

        self.reactant_stoichiometry = {'inhibitor': 0 }

    def modifier(self, expressions):
        """
        change the flux reaction rate expressions
        :param expression: {vnet, vfwd, vbwd}
        :return:
        """
        # Modification of the of Keq
        # expressions = TabDict([('v_net', rate_expression),
        #                        ('v_fwd', forward_rate_expression),
        #                        ('v_bwd', backward_rate_expression),
        #                        ])
        activation = self.get_qssa_rate_expression()
        expressions['v_bwd'] = expressions['v_bwd'] * activation
        expressions['v_fwd'] = expressions['v_fwd'] * activation
        expressions['v_net'] = expressions['v_fwd'] - expressions['v_bwd']

    def get_qssa_rate_expression(self):
        i = self.reactants.inhibitor.symbol
        k = self.parameters.k_inhibition.symbol
        h = self.parameters.hill_coefficient.symbol
        return ( 1 / (1 + (i/k)**h ))

    def update_qssa_rate_expression(self):
        return None

    def get_full_rate_expression(self):
        raise NotImplementedError

    def calculate_rate_constants(self):
        raise NotImplementedError
