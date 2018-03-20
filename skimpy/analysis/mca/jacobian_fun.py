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

class ElasticityFunction:
    def __init__(self, reduced_stoichometry,
                       depednent_elasticity_function,
                       independent_elaticity_function,
                       weights_depednet_metabolites ):

        self.reduced_stoichometry = reduced_stoichometry
        self.depednent_elasticity_function  = depednent_elasticity_function
        self.independent_elaticity_function = independent_elaticity_function
        self.weights_depednet_metabolites = weights_depednet_metabolites


    def __call__(self,fluxes,concentrations,parameters):

        #Calculate the Jacobian
        jacobian = []

        return jacobian
