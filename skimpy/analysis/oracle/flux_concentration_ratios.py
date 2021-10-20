# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2021 Laboratory of Computational Systems Biotechnology (LCSB),
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

from collections import OrderedDict
from copy import deepcopy

from tqdm import tqdm
import pandas as pd
from cobra.util.solver import set_objective
from optlang.exceptions import SolverError

from pytfa.optim.utils import chunk_sum, symbol_sum
from pytfa.optim.variables import ModelVariable, BinaryVariable
from pytfa.optim.constraints import ModelConstraint
from pytfa.utils import numerics

import numpy as np
import sympy

MIN_C = 1e-10
MAX_C = 1e-1

BIGM = numerics.BIGM
EPSILON = numerics.EPSILON

# Namespace
LOG = 'log'
LIN = 'lin'

class BinUseVariable(ModelVariable,BinaryVariable):
    """

    """
    def __init__(self, model, id_, **kwargs):
        BinaryVariable.__init__(self,
                                id_,
                                model=model,
                                **kwargs)
    prefix = 'BinUseVariable_'


class BinVariable(ModelVariable):
    """

    """
    def __init__(self, model, id_, **kwargs):
        ModelVariable.__init__(self,
                                   model,
                                   id_,
                                   **kwargs)
    prefix = 'BinVariable_'


class FluxRatioCons(ModelConstraint):
    """
    Class to represent thermodynamics constraints.
    G: Flux_FW + Fluw_BW > min_flux
    """
    def __init__(self, model, expr, id_, **kwargs):
        ModelConstraint.__init__(self,
                                   model,
                                   expr,
                                   id_,
                                   **kwargs)
    prefix = 'FluxRatioCons_'


def impose_turnover_concentation_ratios(tmodel, metabolites, tva, ratio, in_place=False,
                                    discretization=LOG, N=11):

    if in_place:
        model = tmodel
    else:
        model = tmodel.copy()

    LC_VARS = [l.variable for l in model.log_concentration.get_by_any(metabolites)]

    # for each metabolites get reactions with stoich -1
    mets = model.metabolites.get_by_any(metabolites)

    FLUX_FWD_VARS = []
    FLUX_BWD_VARS = []

    for met in mets:
        this_fdw_flux = 0
        this_bwd_flux = 0
        for rxn in met.reactions:
            if tva.loc[rxn.id,'minimum']*rxn.metabolites[met] > 0 or \
                tva.loc[rxn.id,'maximum']*rxn.metabolites[met] > 0:
                this_fdw_flux += rxn.forward_variable
                this_bwd_flux += rxn.reverse_variable

        FLUX_FWD_VARS.append(this_fdw_flux)
        FLUX_BWD_VARS.append(this_bwd_flux)

    # Add constraints
    for LC, FWD, BWD in zip(LC_VARS, FLUX_FWD_VARS, FLUX_BWD_VARS):
        add_ratio_constraints(model, LC, FWD, BWD, ratio,
                              concentration_range=(np.exp(LC.lb), np.exp(LC.ub)),
                              discretization=discretization,
                              N=N)
        # Check if can be optimized
        print(model.optimize())

    return model



def add_ratio_constraints(model, lc, fwd, bwd, ratio,
                          concentration_range=(MIN_C, MAX_C),
                          discretization=LOG,
                          N=11 ):

    if discretization == LOG:
        bins = np.logspace(np.log10(concentration_range[0]), np.log10(concentration_range[1]), N)
        print(bins)
    elif discretization == LIN:
        bins = np.linspace(concentration_range[0], concentration_range[1], N)
    else:
        raise ValueError('Discretization {} not supported'.format(discretization))


    """
    Constraint the ratio of the flux and the concentration binwise i.e. 
    V/C > ratio -> V/exp(LC) >= ratio  -> V >= ratio*exp(LC)
    
    Approximation:
    V >= ratio* C
    lower_bin <= C <= upper_bin if log(lower_bin) <= LC <= log(upper_bin)
    
    In constraints:
    # For each bin [lower_bin, upper_bin]
    V - ratio* C >= 0 
    C = sum(Ci)
    ub_i * bin_use_i - Ci <= 0
    
    np.log(lb_i) * bin_use_i <= LCi <= np.log(ub_i) * bin_use_i
    LC = sum(LCi)
    
    0 <= sum(bin_use_i) <= 1
    """

    pos_v = fwd + bwd
    C = model.add_variable(BinVariable,
                           id_='concentration_{}'.format(lc.name),
                           hook=model,
                           lb = 0,
                           ub = BIGM,
                           )

    # Ratio constraint
    expr = pos_v - ratio*C
    ratio_cons = model.add_constraint(FluxRatioCons,
                                      id_="ratio_{}".format(lc.name),
                                      hook=model,
                                      expr=expr,
                                      lb=0)
    this_bin_use_variables = []
    this_LCi_variables = []
    this_Ci_variables = []
    for i in range(1, N):
        lower_bin = bins[i-1]
        upper_bin = bins[i]

        bin_use = model.add_variable(BinUseVariable,
                                     id_='ratio_{}_bin_{}_{}'.format(lc.name, lower_bin, upper_bin),
                                     hook=model,)

        this_bin_use_variables.append(bin_use)

        Ci = model.add_variable(BinVariable,
                               id_='concentration_{}_bin_{}_{}'.format(lc.name, lower_bin, upper_bin),
                               hook=model,
                               lb = -BIGM,
                               ub = BIGM,
                               )
        this_Ci_variables.append(Ci)

        # Concentration to log-concentration coupling Ci == UB or Ci == 0
        expr = Ci - upper_bin*bin_use
        bin_coupling = model.add_constraint(FluxRatioCons,
                                            id_='bin_coupling_{}_bin_{}_{}'.format(lc.name, lower_bin, upper_bin),
                                            hook=model,
                                            expr=expr,
                                            lb= -EPSILON,
                                            ub=  EPSILON,)

        LCi = model.add_variable(BinVariable,
                                id_='log_concentration_{}_bin_{}_{}'.format(lc.name, lower_bin, upper_bin),
                                hook=model,
                                lb = -BIGM,
                                ub = BIGM,
                                )
        this_LCi_variables.append(LCi)

        # Bin use
        expr = np.log(lower_bin)*bin_use - LCi
        lower_bin_use = model.add_constraint(FluxRatioCons,
                                       id_='lower_bin_use_{}_bin_{}_{}'.format(LCi.name,lower_bin, upper_bin),
                                       hook=model,
                                       expr=expr,
                                       ub=0,
                                       )
        expr = np.log(upper_bin)*bin_use - LCi
        upper_bin_use = model.add_constraint(FluxRatioCons,
                                   id_='upper_bin_use_{}_bin_{}_{}'.format(LCi.name, lower_bin, upper_bin),
                                   hook=model,
                                   expr=expr,
                                   lb=0,
                                   )

    # Force use only one constraint
    expr = symbol_sum(this_bin_use_variables)
    force_use = model.add_constraint(FluxRatioCons,
                                     id_='force_use_{}'.format(lc.name),
                                     hook=model,
                                     expr=expr,
                                     lb=0,
                                     ub=1,)
    
    # Equality sum(LCi) == LC and sum(Ci) == C
    expr = symbol_sum(this_Ci_variables) - C
    force_use = model.add_constraint(FluxRatioCons,
                                 id_='force_Ci_{}'.format(lc.name),
                                 hook=model,
                                 expr=expr,
                                 lb= -EPSILON,
                                 ub=  EPSILON,)

    expr = symbol_sum(this_LCi_variables) - lc
    force_use = model.add_constraint(FluxRatioCons,
                                 id_='force_LCi_{}'.format(lc.name),
                                 hook=model,
                                 expr=expr,
                                 lb= -EPSILON,
                                 ub=  EPSILON,)
