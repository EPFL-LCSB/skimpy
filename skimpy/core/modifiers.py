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

    def __init__(self, name, modifier = None):
        self._name = name
        if modifier is not None:
            self._modifier = modifier

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

    def __init__(self, name, modifier = None):
        ExpressionModifier.__init__(self, name, modifier)


class ConstantConcentration(BoundaryCondition):
    """

    """

    prefix = 'CC'

    def __init__(self, reactant, name = None):

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

class AdditiveConcentrationRate(ExpressionModifier):
    """
    Add a concentration rate term to your rate expression
    """
    # FIXME Please give us an alternate name we _REALLY_ don't like it

    prefix = 'ADDCR'

    def __init__(self, reactant, flux_value, name=None):

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

    def __init__(self, reactant, flux_value):
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

    def __init__(self, reaction, reactant, stoichiometry, name=None):

        if name is None:
            name = reactant.__repr__()

        reactants = self.Reactants(small_molecule=reactant)
        parameters = self.Parameters()
        KineticMechanism.__init__(self, name, reactants, parameters)

        self.stoichiometry = stoichiometry

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

        if self.stoichiometry < 0:
            expressions['v_fwd'] = expressions['v_fwd']\
                                   * self.get_qssa_rate_expression()

        if self.stoichiometry > 0:
            expressions['v_bwd'] = expressions['v_bwd'] \
                                   * self.get_qssa_rate_expression()

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

    Reactants = make_reactant_set(__name__, ['small_molecule'])

    Parameters = make_parameter_set(    __name__,
                                        { })

    parameter_reactant_links = {}

    def __init__(self, reaction, reactant, stoichiometry, name=None):

        if name is None:
            name = reactant.__str__()

        reactants = self.Reactants(small_molecule=reactant)
        parameters = self.Parameters()
        KineticMechanism.__init__(self, name, reactants, parameters)

        self.stoichiometry = stoichiometry

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
                                * self.get_qssa_rate_expression()**self.stoichiometry

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




