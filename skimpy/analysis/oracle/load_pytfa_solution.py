# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2020 Laboratory of Computational Systems Biotechnology (LCSB),
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

from skimpy.utils.general import sanitize_cobra_vars
from skimpy.utils.conversions import deltag0_to_keq
from skimpy.core.parameters import ParameterValues

import pandas as pd
import numpy as np

# Load and convert pytfa solution for kinetic model

def load_fluxes(solution_raw,tmodel,kmodel,
                density=None,
                ratio_gdw_gww=None,
                concentration_scaling=None,
                time_scaling=None,
                xmol_in_flux=1e-3):
    # TODO try to fetch from model
    if density is None \
        or ratio_gdw_gww is None \
        or concentration_scaling is None \
        or time_scaling is None:
        raise  ValueError("density, ratio_gdw_gww, concentration_scaling, or time_scaling "
                          "is required as input or field of kmodel")

    # Flux solution input assumed to be mmol/gDW/hr
    flux_scaling_factor =  xmol_in_flux * (ratio_gdw_gww * density) \
                           * concentration_scaling \
                           / time_scaling

    fluxes_in_kmodel = list(kmodel.reactions.keys())

    # Convert to net-fluxes
    solution_nf =  { this_rxn.id: (solution_raw[this_rxn.forward_variable.name] \
                      - solution_raw[this_rxn.reverse_variable.name])  \
                     for this_rxn in tmodel.reactions}

    # Convert tmodel net fluxes to kmodel fluxes
    flux_dict = {rxn: solution_nf[rxn]*flux_scaling_factor for rxn in fluxes_in_kmodel}

    fluxes = pd.Series(flux_dict)
    # Sort according to the k-model
    return fluxes[fluxes_in_kmodel]


def load_concentrations(solution_raw, tmodel, kmodel, concentration_scaling=None):
    # TODO try to fetch from model
    if concentration_scaling is None:
        raise  ValueError("concentration_scaling is required as input or field of kmodel")

    concentration_dict = {sanitize_cobra_vars(lc.id): np.exp(solution_raw[lc.variable.name])
                                                      *concentration_scaling
                          for lc in tmodel.log_concentration}
    concentrations = pd.Series(concentration_dict)

    return concentrations


def load_equilibrium_constants(solution_raw, tmodel, kmodel,
                               concentration_scaling=None,
                               in_place=False):
    # TODO try to fetch from model
    if concentration_scaling is None:
        raise  ValueError("concentration_scaling is required as input or field of kmodel")

    equilibrium_constant_dict = dict()

    # Calculate the fitting equilibrium constants for the kinetic models (absorb concentrations that do
    # not appear explicitly in the mass balance of the model in the deltag )
    RT = tmodel.RT

    rnxs_ids_with = [dg.id for dg in tmodel.delta_g]

    for pytfa_rxn in tmodel.reactions:

        if not pytfa_rxn.id in kmodel.reactions:
            continue

        if not pytfa_rxn.id in rnxs_ids_with:
            continue

        deltag0 = solution_raw[tmodel.delta_g.get_by_id(pytfa_rxn.id).name]

        for met, stoich in pytfa_rxn.metabolites.items():
            kin_met_id = sanitize_cobra_vars(met)
            if (kin_met_id in kmodel.reactants) or  (kin_met_id in kmodel.parameters):
                var_met_lc = tmodel.log_concentration.get_by_id(met.id).name
                met_lc = solution_raw[var_met_lc]
                deltag0 -= stoich * RT * (met_lc + np.log(concentration_scaling))

        try:
            k_eq = kmodel.reactions[pytfa_rxn.id].parameters['k_equilibrium']
        except KeyError:
            continue

        if in_place:
            k_eq.value = deltag0_to_keq(deltag0, tmodel.TEMPERATURE,
                                       gas_constant=tmodel.GAS_CONSTANT)
        # Check how to best index this ...
        equilibrium_constant_dict[k_eq.symbol] = deltag0_to_keq(deltag0, tmodel.TEMPERATURE,
                                                                gas_constant=tmodel.GAS_CONSTANT)

    return ParameterValues(equilibrium_constant_dict,kmodel)




