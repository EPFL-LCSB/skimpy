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

from pytfa.optim.variables import ThermoDisplacement, DeltaGstd, DeltaG, \
        ForwardUseVariable, BackwardUseVariable, LogConcentration, \
        ReactionVariable, MetaboliteVariable

from pytfa.optim.constraints import SimultaneousUse, NegativeDeltaG, \
        BackwardDeltaGCoupling, ForwardDeltaGCoupling, BackwardDirectionCoupling, \
        ForwardDirectionCoupling, ReactionConstraint, MetaboliteConstraint, \
        DisplacementCoupling

from math import log

from pytfa.utils import numerics
BIGM = numerics.BIGM
BIGM_THERMO = numerics.BIGM_THERMO
BIGM_DG = numerics.BIGM_DG
BIGM_P = numerics.BIGM_P
EPSILON = numerics.EPSILON
MAX_STOICH = 10


def add_undefined_delta_g(_tmodel,
                          solution,
                          delta_g_std=-10,
                          delta_g_std_err=2,
                          add_displacement=True,
                          inplace=True,
                          exclude_reactions=[]):
    if inplace:
        tmodel = _tmodel
    else:
        tmodel = _tmodel.copy()

    sol = solution

    for this_rxn in tmodel.reactions:
        # Pass over the exclude reactions
        if this_rxn.id in exclude_reactions:
            continue

        this_rxn_rev_flux = sol[this_rxn.id]
        try:
           this_dgo = tmodel.delta_gstd.get_by_id(this_rxn.id)
        except KeyError:
           if not this_rxn.boundary:
                if this_rxn_rev_flux < 0:
                    this_dgo = -delta_g_std
                else:
                    this_dgo = delta_g_std

                add_dummy_delta_g(tmodel,
                                  this_rxn,
                                  delta_g_std=this_dgo,
                                  delta_g_std_err=delta_g_std_err,
                                  add_displacement=add_displacement)
                #Adding dgo
                print(this_rxn.id)
                tmodel.optimize()

    tmodel.repair()

    return tmodel


def add_dummy_delta_g(tmodel,
                      rxn,
                      delta_g_std=-100,
                      delta_g_std_err=2,
                      add_displacement=True):
    RT = tmodel.RT

    DGR_lb = -BIGM_THERMO  # kcal/mol
    DGR_ub = BIGM_THERMO  # kcal/mol


    epsilon = tmodel.solver.configuration.tolerances.feasibility

    """
    Add the delta G as variable
    """
    # add the delta G as a variable
    DGR = tmodel.add_variable(DeltaG, rxn, lb=DGR_lb, ub=DGR_ub)

    # add the delta G naught as a variable
    RxnDGerror = delta_g_std_err
    DGoR = tmodel.add_variable( DeltaGstd,
                                rxn,
                                lb=delta_g_std - RxnDGerror,
                                ub=delta_g_std + RxnDGerror)
    # RxnDGnaught on the right hand side
    RHS_DG = DeltaGstd
    LC_ChemMet = 0

    for met in rxn.metabolites:
        metformula = met.formula
        if metformula not in ['H', 'H2O']:
            # we use the LC here as we already accounted for the
            # changes in deltaGFs in the RHS term
            try:
                tmodel.log_concentration.get_by_id(met.id).variable
            except KeyError:
                metComp = met.compartment
                metLConc_lb = log(tmodel.compartments[metComp]['c_min'])
                metLConc_ub = log(tmodel.compartments[metComp]['c_max'])

                LC = tmodel.add_variable(LogConcentration,
                                       met,
                                       lb=metLConc_lb,
                                       ub=metLConc_ub)
                #tmodel.LC_vars[met] = LC
                tmodel.repair()

            LC_ChemMet += (tmodel.log_concentration.get_by_id(met.id).variable
                           * RT
                           * rxn.metabolites[met])

    # G: - DGR_rxn + DGoRerr_Rxn
    #   + RT * StoichCoefProd1 * LC_prod1
    #   + RT * StoichCoefProd2 * LC_prod2
    #   + RT * StoichCoefSub1 * LC_subs1
    #   + RT * StoichCoefSub2 * LC_subs2
    #   - ...
    #   = 0

    # Formulate the constraint
    CLHS = DGoR - DGR + LC_ChemMet
    tmodel.add_constraint(NegativeDeltaG, rxn, CLHS, lb=0, ub=0)

    if add_displacement:
        lngamma = tmodel.add_variable(ThermoDisplacement,
                                    rxn,
                                    lb=-BIGM_P,
                                    ub=BIGM_P)

        # ln(Gamma) = +DGR/RT (DGR < 0 , rxn is forward, ln(Gamma) < 0d
        expr = lngamma - 1 / RT * DGR
        tmodel.add_constraint(DisplacementCoupling,
                            rxn,
                            expr,
                            lb=0,
                            ub=0)

    # Create the use variables constraints and connect them to the
    # deltaG if the reaction has thermo constraints
    # FU_rxn: 1000 FU_rxn + DGR_rxn < 1000 - epsilon
    FU_rxn = tmodel.add_variable(ForwardUseVariable, rxn)

    CLHS = DGR + FU_rxn * BIGM_THERMO
    tmodel.add_constraint(ForwardDeltaGCoupling,
                        rxn,
                        CLHS,
                        ub=BIGM_THERMO - epsilon*1e3)

    # BU_rxn: 1000 BU_rxn - DGR_rxn < 1000 - epsilon
    BU_rxn = tmodel.add_variable(BackwardUseVariable, rxn)

    CLHS = BU_rxn * BIGM_THERMO - DGR
    tmodel.add_constraint(BackwardDeltaGCoupling,
                        rxn,
                        CLHS,
                        ub=BIGM_THERMO - epsilon*1e3)

    # create the prevent simultaneous use constraints
    # SU_rxn: FU_rxn + BU_rxn <= 1
    CLHS = FU_rxn + BU_rxn
    tmodel.add_constraint(SimultaneousUse, rxn, CLHS, ub=1)

    # create constraints that control fluxes with their use variables
    # UF_rxn: F_rxn - M FU_rxn < 0
    F_rxn = rxn.forward_variable
    CLHS = F_rxn - FU_rxn * BIGM
    tmodel.add_constraint(ForwardDirectionCoupling, rxn, CLHS, ub=0)

    # UR_rxn: R_rxn - M RU_rxn < 0
    R_rxn = rxn.reverse_variable
    CLHS = R_rxn - BU_rxn * BIGM

    tmodel.add_constraint(BackwardDirectionCoupling, rxn, CLHS, ub=0)

