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
from skimpy.core import *

def import_kinetic_mechanism(cobra_reaction,
                             small_molecules = None,
                             water = None ):
    """

    :param cobra_reaction:
    :return: kinetic_mechanism: skimpy kinetic mechanism object
    """
    if small_molecules is None:
        small_molecules = ['h','co2','nh4','o2','pi','ppi','pppi','h2o2','hco3','h2s','so3','so4']
    if water is None:
        water = 'h2o'

    # remove water

    # split in small molecules and reactants
    mets = cobra_reaction.metabolites

    # Get a mechanism
    kinetic_mechanism = None
    # 1) Check for Michaelis-Menten Mechanism

    # 2) Check for Hill Mechanism

    # 3) Convenience kinetics



    # Add small molecules modifiers

    return kinetic_mechanism


