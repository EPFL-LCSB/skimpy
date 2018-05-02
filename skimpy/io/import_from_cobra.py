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
from .import_kinetic_mechanism import import_kinetic_mechanism


def import_from_cobra(cobra_model):
    """
    Function to create a kinetic model from a constraint based model
    :param cobra_model:
    :return: skimpy model
    """

    skimpy_model = KineticModel()

    """Read the reactions """

    # DM_     Boundary reactions

    # By default the metabolites of boundary reactions
    # to be constant concentrations
    for this_reaction in cobra_model.reactions:
        if this_reaction.name.startswith("DM_"):
            for this_met in this_reaction.metabolites:
                this_const_met = ConstantConcentration(this_met)
                skimpy_model.add_boundary_condition(this_const_met)
        else:
            this_kinetic_reaction = import_kinetic_mechanism(this_reaction)
            skimpy_model.add_reaction(this_kinetic_reaction)

    return skimpy_model


