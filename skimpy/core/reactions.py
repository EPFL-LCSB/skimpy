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


class Reaction(object):
    """
    Global reaction class
    """
    def __init__(self, name, substrates, mechanism, parameters=None):
        self.name = name
        self.mechanism = mechanism(name = name,
                                   substrates = substrates,
                                   parameters = parameters)


    # Hooks to the mechanism attributes for convenience
    @property
    def substrates(self):
        return self.mechanism.substrates

    @substrates.setter
    def substrates(self, value):
        self.mechanism.substrates = value

    @property
    def parameters(self):
        return self.mechanism.parameters

    @parameters.setter
    def parameters(self, value):
        self.mechanism.parameters = value

    @property
    def rates(self):
        return self.mechanism.rates

    @rates.setter
    def rates(self, value):
        self.mechanism.rates = value


    def __str__(self):
        return "%s of with %s kinetics "%(self.name, self.type)

    def parametrize(self, params):
        self.mechanism.parameters = params