
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
from cobra.util.solver import set_objective
from pytfa.analysis import variability_analysis
from pytfa.optim.constraints import ReactionConstraint
from pytfa.utils import numerics

class MinFLux(ReactionConstraint):
    """
    Class to represent thermodynamics constraints.
    G: Flux_FW + Fluw_BW > min_flux
    """
    prefix = 'MF_'


BIGM = numerics.BIGM
BIGM_THERMO = numerics.BIGM_THERMO
BIGM_DG = numerics.BIGM_DG
BIGM_P = numerics.BIGM_P
EPSILON = numerics.EPSILON


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
                # TODO if
                expression = fwd_var + rev_var
                temp_model.add_constraints(MinFLux,
                                           rxn,
                                           expression,
                                           lb=flux)
    
    else:
        for rxn in temp_model.reactions:
            if rxn.id not in exclude:
                rev_var = rxn.reverse_variable
                fwd_var = rxn.forward_variable

                expression = fwd_var + rev_var

                temp_model.add_constraints(MinFLux,
                                           rxn,
                                           expression,
                                           lb=flux)

    temp_model.repair()

    return temp_model


