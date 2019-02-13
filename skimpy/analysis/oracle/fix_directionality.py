# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2018 Laboratory of Computational Systems Biotechnology (LCSB),
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

def fix_directionality(tmodel,solution, inplace = True):
    """
    Takes a flux solution and transfers its reaction directionality as
    constraints for the cobra_model
    :param inplace:
    :param tmodel:
    :param solution:
    :return:
    """

    if inplace:
        _tmodel = tmodel
    else:
        _tmodel = tmodel.copy()

    tol = tmodel.solver.configuration.tolerances.feasibility

    for this_reaction in _tmodel.reactions:

        rev_var = this_reaction.reverse_variable
        fwd_var = this_reaction.forward_variable

        if solution.fluxes[this_reaction.id] > -tol:
            rev_var.lb = 0.0
            rev_var.ub = 0.0

        elif solution.fluxes[this_reaction.id] < tol:
            fwd_var.lb = 0.0
            fwd_var.ub = 0.0

    _tmodel.repair()
    return _tmodel
