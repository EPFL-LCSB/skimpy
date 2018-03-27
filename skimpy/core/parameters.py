# -*- coding: utf-8 -*-
"""
.. module:: pytfa
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

from sympy import sympify, Symbol

from ..utils.tabdict import TabDict

class Parameter(object):

    def __init__(self, name, value = None, model=None):

        self.name = name
        self.model = model
        self.value = value
        self._symbol = None
        self._generate_symbol()

    def _generate_symbol(self):
        self._symbol = Symbol(self.name)

    @property
    def symbol(self):
        return self._symbol


class ParameterSet(ABC, TabDict):
    def __init__(self, mechanism, param_names, param_values):

        self.mechanism = mechanism

        TabDict.__init__(self)

        for p,v in param_values.items():
            self[p] = Parameter(p, value = v)

def make_parameter_set(mechanism, param_names):
    name = mechanism.__name__+ParameterSet.__name__

    def this_init(self, **kwargs):
        return ParameterSet.__init__(self,
                                     mechanism=mechanism,
                                     param_names=param_names,
                                     param_values=kwargs)

    return type(name, (ParameterSet,), {'__init__': this_init})