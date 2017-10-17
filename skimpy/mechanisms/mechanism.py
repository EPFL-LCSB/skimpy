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

from collections import namedtuple


class KineticMechanism(ABC):
    def __init__(self, name, substrates, parameters):
        ABC.__init__()
        self.name = name
        self._substrates = substrates
        self._parameters = parameters

        for metafield in [self._substrates, self._parameters]:
            for field in metafield.__dict__.keys():
                setattr(self, field) = property(
                    lambda self:getattr(metafield,field))

    @abstractproperty
    def Substrates(self):
        """
        Class to define metabolites and their roles in the reaction
        :return:
        """

    @abstractproperty
    def Parameters(self):
        """
        Class to define parameters and their roles in the reaction
        :return:
        """
        pass

    @abstractproperty
    def Rates(self):
        """
        Class to define rates and their roles in the reaction
        :return:
        """
        pass

    @abstractmethod
    def calculate_rates(self):
        pass

    @abstractmethod
    def get_qssa_rate_expression(self):
        pass

    @abstractmethod
    def get_full_rate_expression(self):
        pass

