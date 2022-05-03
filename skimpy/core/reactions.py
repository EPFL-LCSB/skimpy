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
from skimpy.utils.namespace import *


class Reaction(object):
    """
    Global reaction class
    """
    def __init__(self, name, reactants, mechanism, parameters=None, inhibitors=None, enzyme=None):
        self.name = name
        self.enzyme = enzyme

        if inhibitors is not None:
            self.mechanism = mechanism(name=name,
                                       reactants=reactants,
                                       parameters=parameters,
                                       inhibitors=inhibitors,
                                       enzyme=enzyme)
        else:
            self.mechanism = mechanism(name=name,
                                       reactants=reactants,
                                       parameters=parameters,
                                       enzyme=enzyme)
        self.modifiers = TabDict([])

    # Hooks to the mechanism attributes for convenience
    @property
    def reactants(self):
        reactants = TabDict( (k,r) for k,r in self.mechanism.reactants.items()
                                    if r.type == VARIABLE)

        for this_modifier in self.modifiers.values():
            this_reactants = TabDict( ("{}_{}".format(k,r.name),r)
                                     for k,r in this_modifier.reactants.items()
                                     if r.type == VARIABLE)
            reactants.update(this_reactants)

        if not self.mechanism.inhibitors is None:
            regulators = TabDict( (k,r) for k,r in self.mechanism.inhibitors.items()
                                 if r.type == VARIABLE)
            reactants.update(regulators)

        return reactants

    @reactants.setter
    def reactants(self, value):
        self.mechanism.reactants = value

    @property
    def reactant_stoichiometry(self):
        reactant_stoichiometry = TabDict([])
        for k,v in self.mechanism.reactant_stoichiometry.items():
            if self.mechanism.reactants[k].type == VARIABLE:
                this_reactant = self.mechanism.reactants[k]
                if this_reactant in reactant_stoichiometry.keys():
                    reactant_stoichiometry[this_reactant] += v
                else:
                    reactant_stoichiometry[this_reactant] = v

        for this_modifier in self.modifiers.values():
            for k,v in this_modifier.reactants.items():
                if v.type == VARIABLE and this_modifier.reactant_stoichiometry[k] != 0:
                    reactant_stoichiometry[v] = this_modifier.reactant_stoichiometry[k]
        return reactant_stoichiometry

    @property
    def parameters(self):
        parameters = TabDict( (k, r) for k, r in self.mechanism.parameters.items()
                                     if r.type == PARAMETER)
        parameters.update(TabDict((k, r) for k, r in self.mechanism.reactants.items()
                             if r.type == PARAMETER))

        for this_modifier in self.modifiers.values():
            # TODO CHECK IF THIS IS A GOOD IDEA
            this_params = TabDict( ( str(r.symbol), r) for k, r in this_modifier.parameters.items()
                                          if r.type == PARAMETER)

            this_params.update(TabDict((k, r) for k, r in this_modifier.reactants.items()
                                              if r.type == PARAMETER))
            parameters.update(this_params)
        return parameters

    @parameters.setter
    def parameters(self, value):
        # Exception for expression based mechanisms
        if not 'ExpressionBasedKinetics' in self.mechanism.__class__.__name__:
            for name, p in value.items():
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