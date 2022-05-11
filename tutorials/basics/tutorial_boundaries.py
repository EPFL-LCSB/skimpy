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

# Test models
from skimpy.core import *
from skimpy.mechanisms import *
from skimpy.utils.namespace import ELEMENTARY


if __name__ == '__main__':

    name = 'pfk'
    metabolites = ReversibleMichaelisMenten.Reactants(substrate='A', product='B')

    ## ELEMENTARY Method
    parameters = ReversibleMichaelisMenten.Parameters(
        vmax_forward=1,
        k_equilibrium=1,
        km_substrate=10,
        km_product=10,
    )

    import warnings
    warnings.warn('We still need to implement vmax to elementary')

    pfk = Reaction(name=name,
                   mechanism=ReversibleMichaelisMenten,
                   reactants=metabolites,
                   )

    this_model = KineticModel()
    this_model.add_reaction(pfk)
    this_model.parametrize_by_reaction({pfk.name:parameters})

    ## Make the Boundary Condition
    the_boundary_condition = ConstantConcentration(this_model.reactants['A'])

    this_model.add_boundary_condition(the_boundary_condition)

    this_model.parameters.A.value = 10.0

    this_model.compile_ode()

    this_model.initial_conditions.B = 1.0

    this_sol_full = this_model.solve_ode([0, 100.0], solver_type='cvode')

    this_sol_full.plot('boundary_out_full.html')
