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

from collections import defaultdict, OrderedDict
import numpy as np

from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.linalg import inv as sparse_inv
from sympy import diff, simplify, Matrix, eye, zeros

from .elasticity_fun import ElasticityFunction

from skimpy.utils.general import get_stoichiometry, join_dicts
from skimpy.utils.tabdict import iterable_to_tabdict, TabDict
from skimpy.utils.namespace import *
from skimpy.utils.moieties import rational_left_basis

from ...utils.namespace import *

sparse_matrix = csc_matrix

def make_mca_functions(kinetic_model,parameter_list,sim_type):
    """ Create the elasticity and flux functions for MCA
    :param kinmodel:
    :param parameter_list:
    :return:
    """

    # Get all variables and expressions (Better solution with types?)
    # TODO This should be a method in KineticModel that stores the expressions
    if sim_type == QSSA:
        all_data = []

        #TODO Modifiers sould be applicable for all simulation types
        for this_reaction in kinetic_model.reactions.values():
            this_reaction.mechanism.get_qssa_rate_expression()
            # Update rate expressions
            for this_mod in this_reaction.modifiers.values():
                this_mod(this_reaction.mechanism.reaction_rates)
            this_reaction.mechanism.update_qssa_rate_expression()

            # Add modifier expressions
            for this_mod in this_reaction.modifiers.values():
                small_mol = this_mod.reactants['small_molecule']
                sm = small_mol.symbol
                flux = this_reaction.mechanism.reaction_rates['v_net']
                flux_expression_sm = flux*this_mod.stoichiometry
                this_reaction.mechanism.expressions[sm] = flux_expression_sm
                # Add small molecule parameters if they are
                if small_mol.type == PARAMETER:
                    this_reaction.mechanism.expression_parameters.update([small_mol.name])

            all_data.append((this_reaction.mechanism.expressions,
                             this_reaction.mechanism.expression_parameters))

    elif sim_type == TQSSA:
        raise(NotImplementedError)

    elif sim_type == ELEMENTARY:
        raise(NotImplementedError)

    else:
        raise ArgumentError('{} is not recognized as a simulation type'.format(sim_type))


    # Get flux expressions for the net
    all_flux_expressions = [this_reaction.mechanism.reaction_rates['v_net'] \
                           for this_reaction in kinetic_model.reactions.values()]

    all_expr, all_parameters = list(zip(*all_data))

    # Flatten all the lists
    flatten_list = lambda this_list: [item for sublist in this_list \
                                      for item in sublist]

    all_rates = flatten_list([these_expressions.keys()
                              for these_expressions in all_expr])

    # Sort into an ordered list
    all_parameters = flatten_list(all_parameters)
    all_parameters = list(set(all_parameters))
    all_parameters = iterable_to_tabdict(all_parameters, use_name=False)


    # # Get unique set of all the variables
    # all_variables = set(all_rates)
    # all_variables = iterable_to_tabdict(all_variables, use_name=False)

    # Better since this is implemented now
    all_variables = TabDict([(k,v.symbol) for k,v in kinetic_model.reactants.items()])

    #Get depedent variables (i.e. concentrations)
    reduced_stoichiometry, dependent_weights, independent_ix, dependent_ix = \
        get_reduced_stoichiometry(kinetic_model, all_variables)

    all_independent_variables = OrderedDict([(k,v)
                                 for e,(k,v) in enumerate(all_variables.items())
                                 if e in independent_ix])

    all_dependent_variables = OrderedDict([(k,v)
                                 for e,(k,v) in enumerate(all_variables.items())
                                 if e in dependent_ix])


    #parameter elasticity function
    if parameter_list:
        parameter_elasticities_fun = make_elasticity_fun(all_flux_expressions,
                                                         parameter_list,
                                                         all_variables,
                                                         all_parameters)
    else:
        parameter_elasticities_fun = None

    #concentration elasticity functions
    independent_elasticity_fun = make_elasticity_fun(all_flux_expressions,
                                                      all_independent_variables,
                                                      all_variables,
                                                      all_parameters)

    if all_dependent_variables:
        dependent_elasticity_fun = make_elasticity_fun(all_flux_expressions,
                                                        all_dependent_variables,
                                                        all_variables,
                                                        all_parameters)
    else:
        dependent_elasticity_fun = None


    return reduced_stoichiometry,\
           independent_elasticity_fun, \
           dependent_elasticity_fun, \
           parameter_elasticities_fun, \
           dependent_weights,\
           all_variables, \
           independent_ix, \
           dependent_ix




def make_elasticity_fun(expressions,respective_variables ,variables, parameters):
    """
    Create an ElasticityFunction with elasticity = dlog(expression)/dlog(respective_variable)
    :param expressions  tab_dict of expressions (e.g. forward and backward fluxes)
    :param variables    list of variables as string (e.g. concentrations or parameters)
    
    """

    # Get the derivative of expression x vs variable y
    elasticity_expressions = {}
    column = 0
    row = 0

    for this_expression in expressions:

        column = 0
        for this_variable in respective_variables.values():

            this_elasticity = get_dlogx_dlogy(this_expression, this_variable)

            if this_elasticity != 0:
                elasticity_expressions[(row, column)] = this_elasticity
            column += 1
        row += 1

    # Shape of the matrix
    shape = (len(expressions), len(respective_variables))

    # Create the elasticity function
    elasticity_fun = ElasticityFunction(elasticity_expressions, variables, parameters, shape)

    return elasticity_fun


def get_dlogx_dlogy(sympy_expression, variable):
    """
    Calc d_log_x/d_log_y = y/x*dx/dy
    """
    partial_derivative = diff(sympy_expression, variable)

    expression = simplify(partial_derivative / sympy_expression * variable)

    return expression


def get_reduced_stoichiometry(kinetic_model, all_variables):

    #TODO IMPLEMENT THE MOIJETY DETECTION
    full_stoichiometry = get_stoichiometry(kinetic_model, all_variables)

    S = Matrix(full_stoichiometry.todense())
    # Left basis dimensions: rows are metabolites, columns are moieties
    # L0*S = 0 -> L0 is the left null space matrix
    left_basis = rational_left_basis(S)
    L0 = Matrix([x.transpose() for x in left_basis])

    ## We need to separate N and N0 beforehand

    # Per moiety, select one variable that has not been selected before
    # L0
    L0_sparse = sparse_matrix(np.array(L0), dtype=np.float)

    nonzero_rows, nonzero_cols = L0_sparse.nonzero()
    row_dict = defaultdict(list)

    # Put the ixs in a dict indexed by row number (moiety index)
    for k, v in zip(nonzero_rows, nonzero_cols):
        row_dict[k].append(v)

    # The first independent variables are those involved in no moieties
    all_independent_ix = [x for x in range(L0.shape[1]) if not x in nonzero_cols]
    # Indices for dependent metabolites indices
    all_dependent_ix = []

    # For each line, get an exclusive representative.
    # There should be at least as many exclusive representatives as lines

    # Iterate over mojeties and start with the ones with least members
    for row in sorted(row_dict, key=lambda k: len(row_dict[k])):
        mojetie_vars = row_dict[row]
        # Get all unassigned metabolites participating in this mojetie
        unassigned_vars = [x for x in set(mojetie_vars)
            .difference(all_independent_ix+all_dependent_ix)]
        # Choose a representative dependent metabolite:
        if unassigned_vars:
            all_dependent_ix.append(unassigned_vars[-1])
        else:
            raise Exception('Could not find an dependent var that is not already used'
                            ' in {}'.format(mojetie_vars))

    # The independent mets is the set difference from the dependent
    all_independent_ix = [x for x in set(range(L0.shape[1]))
        .difference(all_dependent_ix)]

    # Reindex S in N, N0
    S = S[all_independent_ix+all_dependent_ix,:]
    # If we reindex S, then so should be L0
    L0 = L0[:,all_independent_ix+all_dependent_ix]

    # Getting the reduced Stoichiometry:
    # S is the full stoichiometric matrix
    # N  is the full rank reduced stoichiometric matrix
    # N0 is the remainder
    #
    #     [ I_n (nxn) | 0_r nxr) ]
    # L = [     L0 ((r)x(n+r))    ]
    #
    #     [ N  ]
    # S = [ N0 ]
    #
    # Then:
    #          [ I_n (nxn) | 0_r (nxr) ] * [ N  ]
    # L * S  = [     L0 ((r)x(n+r))    ]   [ N0 ]
    #
    #          [ I_n*N + 0_r * N0 ]   [ N ]
    # L * S  = [      L0 * S      ] = [ 0 ]

    r,n_plus_r = L0.shape
    n = n_plus_r - r

    # Upper block
    U = Matrix([eye(n), zeros(r,n)]).transpose()
    L = Matrix([U,L0])

    N_ = L*S
    N  = N_[:n,:] # the rows after n are 0s

    reduced_stoichiometry   = sparse_matrix(np.array(N),dtype=np.float)
    conservation_relation = L0_sparse

    # If the left hand null space is empty no mojeties
    if not any(L0_sparse):
        reduced_stoichiometry = full_stoichiometry
        all_independent_ix = range(full_stoichiometry.shape[0])
        all_dependent_ix = []

    return reduced_stoichiometry, conservation_relation, all_independent_ix, all_dependent_ix




