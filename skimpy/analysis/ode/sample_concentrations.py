import numpy as np
from cobra import Model, Reaction, Metabolite
#from optlang import Variable, Constraint, Objective, Model
from cobra.sampling import sample
from scipy.sparse import diags

EPSILON = 1e-7

def sample_initial_concentrations(kmodel,
                                  reference_concentrations,
                                  lower_bound=0.8,
                                  upper_bound=1.2,
                                  n_samples=10,
                                  absolute_bounds = False):

    concentrations = np.array([reference_concentrations[k] for k in kmodel.variables])

    if kmodel.conservation_relation is None:
        N = len(kmodel.variables)
        rand = np.random.uniform(low=lower_bound,
                                 high=upper_bound,
                                 size=(N,n_samples))

    #
    else:
        # Todo Also allow
        lower_bound = {k: v * lower_bound for k,v in  reference_concentrations.items()}
        upper_bound = {k: v * upper_bound for k, v in reference_concentrations.items()}


        # Account for volume differences in compartments
        if kmodel.volume_ratio_func is None:
            effective_conservation_relations = kmodel.conservation_relation.todense()
        else:
            param_values = {p.symbol: p.value for p in kmodel.parameters.values()}
            volume_ratios = kmodel.volume_ratio_func(param_values)
            inv_volume_ratio_matrix =  diags(1./np.array(volume_ratios)).todense()
            effective_conservation_relations = kmodel.conservation_relation.todense()\
                                               .dot(inv_volume_ratio_matrix)

        # Create a linear problem and sample it
        rhs = effective_conservation_relations.dot(concentrations)
        linmodel = create_linear_model(effective_conservation_relations,
                                       rhs,
                                       kmodel.reactants,
                                       lower_bound=lower_bound,
                                       upper_bound=upper_bound)

        return sample(linmodel, n_samples, )


def create_linear_model(A, rhs, variables, lower_bound=None, upper_bound=None):
    # This is a bit retarded but this way we can use the sampling functions from cobra

    model = Model('lin_model')
    # TODO have cases for lower and upper bound  None
    rxns = np.array([ Reaction(v, lower_bound=lower_bound[v],
                                  upper_bound=upper_bound[v]) for v in variables])

    model.add_reactions(rxns)

    lin_variables = np.array([r.forward_variable for r in rxns])
    NR, NC = A.shape
    constraints = [model.problem.Constraint( (A[i,:].dot(lin_variables))[0,0] ,
                                             lb=np.float(rhs[:,i]),
                                             ub=np.float(rhs[:,i])) for i in range(NR)]
    model.add_cons_vars(constraints)

    return model
