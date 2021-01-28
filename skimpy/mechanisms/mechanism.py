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
from abc import ABC, abstractmethod
from skimpy.core.itemsets import Reactant
from skimpy.utils.namespace import *


class KineticMechanism(ABC):

    parameter_reactant_links = {}

    def __init__(self,
                 name,
                 reactants,
                 parameters=None,
                 inhibitors=None,
                 enzyme=None):

        # ABC.__init__()
        self.name = name
        self.reactants = reactants
        self.inhibitors = None
        self._parameters = None
        self.enzyme = enzyme

        if parameters is not None:
            self._parameters = parameters
            # self.set_dynamic_attribute_links(self._parameters)
        if inhibitors is not None:
            self.inhibitors = inhibitors

        if enzyme is not None:
            self.reactants['enzyme'] = Reactant(enzyme)

    def __reduce__(self):
            return KineticMechanism.__class__.__name__

    def link_parameters_and_reactants(self):
        for p,r in self.parameter_reactant_links.items():
            try:
                reactant = self.reactants[r]
            except KeyError:
                reactant = self.inhibitors[r]
            parameter = self.parameters[p]
            parameter.hook = reactant
            #reactant.hook = parameter


    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, value):
        self._parameters = value
        self.link_parameters_and_reactants()

    @property
    @abstractmethod
    def Reactants(self):
        """
        Class to define metabolites and their roles in the reaction
        :return:
        """

    @property
    @abstractmethod
    def Parameters(self):
        """
        Class to define parameters and their roles in the reaction
        :return:
        """
        pass


    @abstractmethod
    def get_qssa_rate_expression(self):
        pass

    @abstractmethod
    def update_qssa_rate_expression(self):
        pass

    @abstractmethod
    def get_full_rate_expression(self):
        pass

    @abstractmethod
    def calculate_rate_constants(self):
        pass

    def get_parameters_from_expression(self, expr):

        reactants = [x.symbol for x in self.reactants.values()
                               if x.type == VARIABLE]
        if self.inhibitors is not None:
            inhibitors = [x.symbol for x in self.inhibitors.values()
                         if x.type == VARIABLE]
            reactants += inhibitors

        parameters = set(expr.free_symbols).difference(reactants)

        return parameters


class ElementrayReactionStep(object):
    def __init__(self,educts,products,rate_constant_name):
        self.educts = educts
        self.products = products
        self.rate_constant_name = rate_constant_name

    def __str__(self):
        is_first = True
        for educt in self.educts:
            if is_first:
                educt_string = educt
            else:
                educt_string += " + "+educt
            is_first = False

        is_first = True
        for product in self.products:
            if is_first:
                product_string = product
            else:
                product_string += " + "+product
            is_first = False

        return educt_string + " --> " + product_string

    def __repr__(self):
        return self.__str__()
