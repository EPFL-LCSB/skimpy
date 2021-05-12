
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


from collections import OrderedDict
from copy import deepcopy

from tqdm import tqdm
import pandas as pd
from cobra.util.solver import set_objective
from optlang.exceptions import SolverError

from pytfa.optim.config import dg_relax_config
from pytfa.optim.utils import get_solution_value_for_variables, chunk_sum, symbol_sum
from pytfa.optim.variables import PosSlackVariable, NegSlackVariable, ReactionVariable
from pytfa.utils import numerics

BIGM = numerics.BIGM
BIGM_THERMO = numerics.BIGM_THERMO
BIGM_DG = numerics.BIGM_DG
BIGM_P = numerics.BIGM_P
EPSILON = numerics.EPSILON

class MinFLuxVariable(ReactionVariable):
    """
    Class to represent a negative slack variable for relaxation problems
    """
    def __init__(self,reaction,**kwargs):
        ReactionVariable.__init__(self, reaction, **kwargs)
    prefix = 'MinFluxVar_'

class MinFLux(ReactionConstraint):
    """
    Class to represent thermodynamics constraints.
    G: Flux_FW + Fluw_BW > min_flux
    """
    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)
    prefix = 'MF_'


BIGM = numerics.BIGM
BIGM_THERMO = numerics.BIGM_THERMO
BIGM_DG = numerics.BIGM_DG
BIGM_P = numerics.BIGM_P
EPSILON = numerics.EPSILON


def add_min_flux_requirements(tmodel, flux, inplace=True, exclude=[]):
    if inplace:
        temp_model = tmodel
    else:
        temp_model = deepcopy(tmodel)
        temp_model.repair()

    for rxn in temp_model.reactions:
        # Check input
        if isinstance(flux,dict):
            this_flux = flux[rxn.id]
        else:
            this_flux = flux

        if rxn.id not in exclude:
            rev_var = rxn.reverse_variable
            fwd_var = rxn.forward_variable

            mf_var = temp_model.add_variable(MinFLuxVariable,
                                             rxn,
                                             lb=this_flux,
                                             ub=BIGM,
                                             queue=False)

            expression = fwd_var + rev_var - mf_var
            temp_model.add_constraint( MinFLux,
                                       rxn,
                                       expression,
                                       lb=0,
                                       ub=BIGM,
                                       queue=False)

    temp_model.repair()

    return temp_model



def relax_min_flux(tmodel, reactions_to_ignore=(), solver=None, in_place = False):
    """
    :param t_tmodel:
    :type t_tmodel: pytfa.thermo.ThermoModel:
    :param reactions_to_ignore: Iterable of reactions that should not be relaxed
    :param solver: solver to use (e.g. 'optlang-glpk', 'optlang-cplex',
        'optlang-gurobi'
    :return: a cobra_model with relaxed bounds on standard Gibbs free energy
    """

    if solver is None:
        solver = tmodel.solver.interface

    # Create a copy of the cobra_model on which we will perform the slack addition
    slack_model = deepcopy(tmodel)
    slack_model.solver = solver
    slack_model.name = 'SlackModel '+tmodel.name
    slack_model.id = 'SlackModel_'+tmodel.id

    # Ensure the lazy updates are all done
    slack_model.repair()

    if not in_place:
        # Create a copy that will receive the relaxation
        relaxed_model = deepcopy(tmodel)
        relaxed_model.solver = solver
        relaxed_model.name = 'RelaxedModel '+tmodel.name
        relaxed_model.id = 'RelaxedModel_'+tmodel.id
        relaxed_model.repair()
    else:
        relaxed_model = slack_model

    dg_relax_config(slack_model)

    original_objective = relaxed_model.objective


    # Do not relax if cobra_model is already optimal
    # try:
    #     solution = tmodel.optimize()
    # except SolverError as SE:
    #     status = tmodel.solver.status
    #     tmodel.logger.error(SE)
    #     tmodel.logger.warning('Solver status: {}'.format(status))
    # if tmodel.solver.status == OPTIMAL:
    #     raise Exception('Model is already optimal')

    # Find variables that represent standard Gibbs Energy change
    my_dgo = relaxed_model.get_variables_of_type(MinFLuxVariable)

    # Find constraints that represent negativity of Gibbs Energy change
    my_neg_dg = slack_model.get_constraints_of_type(MinFLux)

    changes = OrderedDict()
    objective_symbols = []

    slack_model.logger.info('Adding slack constraints')

    slack_model.solver.update()

    for this_neg_dg in tqdm(my_neg_dg, desc='adding slacks'):

        # If there is no thermo, or relaxation forbidden, pass
        if this_neg_dg.id in reactions_to_ignore or this_neg_dg.id not in my_dgo:
            continue

        # Create the negative and positive slack variables
        # We can't queue them because they will be in an expression to declare
        # the constraint
        neg_slack = slack_model.add_variable(NegSlackVariable,
                                             this_neg_dg.reaction, lb=0,
                                             ub=BIGM_DG,
                                             queue=False)
        pos_slack = slack_model.add_variable(PosSlackVariable,
                                             this_neg_dg.reaction, lb=0,
                                             ub=BIGM_DG,
                                             queue=False)

        # Create the new constraint by adding the slack variables to the
        # negative delta G constraint (from the initial cobra_model)
        new_expr = this_neg_dg.constraint.expression
        new_expr -= (pos_slack - neg_slack)

        this_reaction = this_neg_dg.reaction
        # # Remove former constraint to override it
        # slack_model.remove_constraint(slack_model._cons_dict[this_neg_dg.name])
        #
        # # Add the new variant
        # slack_model.add_constraint(NegativeDeltaG,
        #                            this_reaction,
        #                            expr=new_expr,
        #                            lb=0,
        #                            ub=0,
        #                            queue=True)

        this_neg_dg.change_expr(new_expr)

        # Update the objective with the new variables
        objective_symbols += [neg_slack,  pos_slack]

    # objective = chunk_sum(objective_symbols)
    objective = symbol_sum(objective_symbols)

    # Change the objective to minimize slack
    set_objective(slack_model, objective)

    # Update variables and constraints references
    slack_model.repair()

    slack_model.logger.info('Optimizing slack model')
    # Relax
    slack_model.objective.direction = 'min'
    relaxation = slack_model.optimize()

    # Extract the relaxation values from the solution, by type
    relaxed_model.logger.info('Extracting relaxation')
    my_neg_slacks = slack_model.get_variables_of_type(NegSlackVariable)
    my_pos_slacks = slack_model.get_variables_of_type(PosSlackVariable)

    neg_slack_values = get_solution_value_for_variables(relaxation,
                                                        my_neg_slacks)
    pos_slack_values = get_solution_value_for_variables(relaxation,
                                                        my_pos_slacks)

    epsilon = relaxed_model.solver.configuration.tolerances.feasibility
    relaxed_model.repair()
    relaxed_model.solver.update()
    # Apply reaction delta G standard bound change
    for this_reaction in tqdm(relaxed_model.reactions, desc = 'applying slack'):
        # No thermo, or relaxation forbidden
        if this_reaction.id in reactions_to_ignore or this_reaction.id not in my_dgo:
            continue

        # Get the standard delta G variable
        the_dgo = my_dgo.get_by_id(this_reaction.id)

        # Get the relaxation
        dgo_delta_lb = \
            neg_slack_values[my_neg_slacks \
                .get_by_id(this_reaction.id).name]
        dgo_delta_ub = \
            pos_slack_values[my_pos_slacks \
                .get_by_id(this_reaction.id).name]

        if in_place:
            the_neg_slack = my_neg_slacks.get_by_id(this_reaction.id)
            the_neg_slack_value = slack_model.solution.raw[the_neg_slack.name]
            the_neg_slack.variable.lb = the_neg_slack_value - epsilon
            the_neg_slack.variable.ub = the_neg_slack_value + epsilon

            the_pos_slack = my_pos_slacks.get_by_id(this_reaction.id)
            the_pos_slack_value = slack_model.solution.raw[the_pos_slack.name]
            the_pos_slack.variable.lb = the_pos_slack_value - epsilon
            the_pos_slack.variable.ub = the_pos_slack_value + epsilon

        # Apply reaction delta G standard bound change
        if dgo_delta_lb > 0 or dgo_delta_ub > 0:

            # Store previous values
            previous_dgo_lb = the_dgo.variable.lb
            previous_dgo_ub = the_dgo.variable.ub

            if not in_place:
                # Apply change
                the_dgo.variable.lb -= (dgo_delta_lb + epsilon)
                the_dgo.variable.ub += (dgo_delta_ub + epsilon)

            # If needed, store that in a report table
            changes[this_reaction.id] = [
                previous_dgo_lb,
                previous_dgo_ub,
                dgo_delta_lb,
                dgo_delta_ub,
                the_dgo.variable.lb,
                the_dgo.variable.ub]


    relaxed_model.repair()
    relaxed_model.logger.info('Testing relaxation')

    relaxed_model.objective = original_objective
    relaxed_model.objective.direction = 'max'

    relaxed_model.optimize()

    if len(changes) == 0:
        # The model is infeasible or something went wrong
        tmodel.logger.error('Relaxation could not complete '
                            '(no DeltaG relaxation found)')
        return relaxed_model, slack_model, None

    # Format relaxation
    relax_table = pd.DataFrame.from_dict(changes,
                                         orient = 'index')
    relax_table.columns = [ 'lb_in',
                            'ub_in',
                            'lb_change',
                            'ub_change',
                            'lb_out',
                            'ub_out']

    return relaxed_model, slack_model, relax_table
