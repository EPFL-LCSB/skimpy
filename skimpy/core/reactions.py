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

from skimpy.utils.tabdict import TabDict

class Reaction(object):
    """
    Global reaction class
    """
    def __init__(self, name, reactants, mechanism, parameters=None):
        self.name = name
        self.mechanism = mechanism(name = name, reactants = reactants,
                                   parameters = parameters)
        self.modifiers = TabDict([])


    # Hooks to the mechanism attributes for convenience
    @property
    def reactants(self):
        reactants = self.mechanism.reactants
        for this_modifier in self.modifiers.values():
            reactants.update(this_modifier.reactants)
        return reactants

    @reactants.setter
    def reactants(self, value):
        self.mechanism.reactants = value

    @property
    def reactant_stoichiometry(self):
        reactant_stoichiometry = TabDict( (self.mechanism.reactants[k],v)
                            for k,v in self.mechanism.reactant_stoichiometry.items())

        for this_modifier in self.modifiers.values():
            this_mod_name = this_modifier.reactants.small_molecule
            reactant_stoichiometry[this_mod_name] = this_modifier.stoichiometry
        return reactant_stoichiometry

    @property
    def parameters(self):
        parameters = self.mechanism.parameters
        for this_modifier in self.modifiers.values():
            parameters.update(this_modifier.parameters)
        return parameters

    # TODO implement the setter
    @parameters.setter
    def parameters(self, value):
        for name,p in value.items():
            p.suffix = self.name
        self.mechanism.parameters = value

    @property
    def rates(self):
        return self.mechanism.rates

    @rates.setter
    def rates(self, value):
        self.mechanism.rates = value


    def __str__(self):
        return "%s of with %s kinetics "%(self.name, self.mechanism)

    def parametrize(self, params):
        self.parameters = params