
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

from pytfa.analysis import variability_analysis


def add_min_flux_requirements(tmodel,flux, inplace=True, safe=True, exclude=[]):
    if inplace:
        temp_model = tmodel
    else:
        temp_model = tmodel.copy()
        temp_model.repair()

    if safe:
        # Make variability analysis to check which reaction can carry the flux
        tva_fluxes = variability_analysis(temp_model, kind='reactions')

        blocked_rxns = tva_fluxes[(tva_fluxes['minimum'] > -flux ) &
                                  (tva_fluxes['maximum'] < flux  ) ].index

        blocked_rxns = [b for b in blocked_rxns if b not in exclude]
        # Remove the reactions that cant carry the minimum flux requirement
        temp_model.remove_reactions([temp_model.reactions.get_by_id(rxn)
                                 for rxn in blocked_rxns])

        for rxn in temp_model.reactions:
            if rxn.id not in exclude:
                rev_var = rxn.reverse_variable
                fwd_var = rxn.forward_variable

                if flux <= rev_var.ub \
                    and rev_var.lb < flux \
                    and tva_fluxes['minimum'][rxn.id] < -flux:
                    rev_var.lb = flux
                if flux <= fwd_var.ub \
                    and fwd_var.lb < flux \
                    and tva_fluxes['maximum'][rxn.id] > flux:
                    fwd_var.lb = flux
    
    else:
        for rxn in temp_model.reactions:
            if rxn.id not in exclude:
                rev_var = rxn.reverse_variable
                fwd_var = rxn.forward_variable

                if flux <= rev_var.ub and rev_var.lb < flux:
                    rev_var.lb = flux
                if flux <= fwd_var.ub and fwd_var.lb < flux:
                    fwd_var.lb = flux

    temp_model.repair()

    return temp_model

def relax_min_flux_gurobi(model, relax_obj_type = 0, exclude_variables=None):
    N = len(model.reactions)
    # Can we exclude the biomass reaction?
    if exclude_variables is None:
        the_vars = [x._internal_variable for x in model.variables[:2*N]]
    else:
        the_vars = [x._internal_variable for x in model.variables[:2 * N]
                    if x not in exclude_variables ]

    vars_penalities = [1]*len(the_vars)

    grm = model.solver.problem.feasRelax(relaxobjtype=relax_obj_type,
                                         minrelax=True,
                                         constrs=None,
                                         vars=the_vars,
                                         lbpen=vars_penalities,
                                         ubpen=None,
                                         rhspen=None)

    return grm

# TODO Flux relaxation function for any solver see pytfa optim relax_LC!

