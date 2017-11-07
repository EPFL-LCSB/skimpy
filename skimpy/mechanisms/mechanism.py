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

class KineticMechanism(ABC):
    def __init__(self, name, substrates, parameters = None):
        # ABC.__init__()
        self.name = name
        self.substrates = substrates
        # self.set_dynamic_attribute_links(self._substrates)

        if parameters is not None:
            self.parameters = parameters
            # self.set_dynamic_attribute_links(self._parameters)

    #
    # def set_dynamic_attribute_links(self, metafield):
    #     for field in metafield.__dict__.keys():
    #         setattr(self, field) = property(
    #             lambda self: getattr(metafield, field))


    @property
    @abstractmethod
    def Substrates(self):
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

    @property
    @abstractmethod
    def RateConstants(self):
        """
        Class to define rates and their roles in the reaction
        :return:
        """
        pass

    @abstractmethod
    def get_qssa_rate_expression(self):
        pass

    @abstractmethod
    def get_full_rate_expression(self):
        pass

    @abstractmethod
    def calculate_rate_constants(self):
        pass

    def get_expression_parameters_from(self,kind):
        parameters = {}
        values_from_kind = getattr(self,kind)

        for the_param in values_from_kind._fields:
            param_value = getattr(values_from_kind,the_param)
            if param_value is None:
                continue
            param_name = the_param + '_' + self.name
            parameters[param_name] = param_value

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
