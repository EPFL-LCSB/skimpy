# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2019 Laboratory of Computational Systems Biotechnology (LCSB),
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

from skimpy.utils.tabdict import TabDict
from sympy import Symbol
from pandas import Series

class ParameterValues(object):
    """
    Parameters set for kinetic models wichan be indexed with symbols or
    """
    def __init__(self,kmodel,parameter_values):
        """

        :param kmodel: KineticModel class
        :param parameter_values: a dict contaning parameter names and values
        """
        if parameter_values.__class__ is Series:
            parameter_values = parameter_values.to_dict()

        self._parameter_values = TabDict([(str(p),v) for p,v in parameter_values.items()])
        # Check if this is a good solution
        model_params = kmodel.parameters

        self._sym_to_str = {model_params[p].symbol:p for p in self._parameter_values}
        self._str_to_param = {p:model_params[p] for p in self._parameter_values}

    def __getitem__(self, item):
        if item.__class__ is  str:
            return self._parameter_values[item]
        if item.__class__ is Symbol:
            return self._parameter_values[self._sym_to_str[item]]

    def items(self):
        return self._parameter_values.items()

    def keys(self):
        return  self._parameter_values.keys()

    def values(self):
        return  self._parameter_values.values()