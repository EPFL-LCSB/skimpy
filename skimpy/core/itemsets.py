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
from ..utils.namespace import *
from warnings import warn

"""

These classes are bases on which each mechanism will build their own 
ParameterSet and ReactantSet, using a trick explained in
https://stackoverflow.com/questions/18402975/python-programatically-create-subclasses-based-on-init-values

The custom ParameterSet/SubstrateSet constructor is called with a custom init 
made on the fly when calling make_parameter_set or make_substrate_set. The 
keyword arguments define the items (Parameters or Substrates) that will be in 
the set. The make_* functions return the custom classes constructors 

"""


class Item(object):

    def __init__(self, name, value=None, model=None, suffix=''):

        self.name = name
        self.model = model
        self.value = value
        self._suffix = suffix
        self._symbol = None
        self._generate_symbol()
        self.hook = None
        self.type = None
        self._lower_bound = None
        self._upper_bound = None

    def _generate_symbol(self):
        if self.suffix:
            self._symbol = Symbol('{}_{}'.format(self.name, self.suffix))
        else:
            self._symbol = Symbol(self.name)

    @property
    def symbol(self):
        return self._symbol

    @property
    def suffix(self):
        return self._suffix

    @suffix.setter
    def suffix(self, value):
        self._suffix = value
        self._generate_symbol()

    @property
    def bounds(self):
        return (self._lower_bound, self._upper_bound)

    @bounds.setter
    def bounds(self, value):
        lower_bound, upper_bound = value
        self._lower_bound = lower_bound
        self._upper_bound = upper_bound


    def __str__(self):
        return str(self.symbol)

    def __repr__(self):
        repr = "{} object {} with symbol {}".format(self.__class__.__name__,
                                                    self.name,
                                                    self.symbol)
        return repr


class ItemSet(ABC, TabDict):
    def __init__(self, mechanism):

        self.mechanism = mechanism

        TabDict.__init__(self)

    def __reduce__(self):
        return self.__class__.__name__ + self.mechanism.__class__.__name__
        

class Parameter(Item):
    def __init__(self, name, required_for=None, value=None, model=None, suffix=''):
        Item.__init__(self, name, value=value, model=model, suffix=suffix)
        self.type = PARAMETER

        if required_for is not None:
            self._required_for = set(required_for)
        else:
            self._required_for = set()


class ParameterSet(ItemSet):
    def __init__(self, mechanism, param_declaration, param_values, suffix=''):

        ItemSet.__init__(self, mechanism=mechanism)

        for p,req in param_declaration.items():
            if p in param_values:
                value = param_values[p]
            else:
                value=None
            self[p] = Parameter(p, required_for=req, value=value, suffix=suffix)

    @property
    def required_for(self):
        return self._required_for


def make_parameter_set(mechanism, param_declaration):
    name = mechanism + ParameterSet.__name__

    def this_init(self, **kwargs):
        return ParameterSet.__init__(self,
                                     mechanism=mechanism,
                                     param_declaration=param_declaration,
                                     param_values=kwargs)

    return type(name, (ParameterSet,), {'__init__': this_init})


class Reactant(Item):
    def __init__(self, name, value = None, model=None, suffix = ''):
        Item.__init__(self, name, value=value, model=model, suffix=suffix)
        self.type = VARIABLE
        self.compartment = None


class ReactantSet(ItemSet):
    def __init__(self, mechanism, reactant_declaration, reactant_values):

        ItemSet.__init__(self, mechanism=mechanism)

        if set(reactant_declaration) != set(reactant_values.keys()):
            if set(reactant_values.keys()).issuperset(reactant_declaration):
                warn('More reactants provided than declared in the signature of '
                               'the mechanism - {} ignored'
                                .format(set(reactant_values).difference(reactant_declaration)))
            else:
                raise KeyError('The reactants provided do not match the signature of '
                               'the mechanism - should contain: {} not {}'
                                .format(set(reactant_declaration),
                                        set(reactant_values)))

        for s,v in reactant_values.items():
            self[s] = Reactant(v)
            # TODO should this be in the Reactant init? or in ItemsSet
            self[s].type = VARIABLE


def make_reactant_set(mechanism, reactant_declaration):
    name = mechanism + ReactantSet.__name__

    def this_init(self, **kwargs):
        return ReactantSet.__init__(self,
                                     mechanism=mechanism, 
                                     reactant_declaration=reactant_declaration,
                                     reactant_values=kwargs)

    return type(name, (ReactantSet,), {'__init__': this_init})
