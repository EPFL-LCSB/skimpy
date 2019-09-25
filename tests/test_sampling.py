import pytest
# Test models
from skimpy.core import *
from skimpy.mechanisms import *
from skimpy.core.solution import ODESolutionPopulation
from skimpy.sampling.simple_parameter_sampler import SimpleParameterSampler


def build_linear_pathway_model():
    # Build linear Pathway model
    metabolites_1 = ReversibleMichaelisMenten.Reactants(substrate='A', product='B')
    metabolites_2 = ReversibleMichaelisMenten.Reactants(substrate='B', product='C')
    metabolites_3 = ReversibleMichaelisMenten.Reactants(substrate='C', product='D')

    ## QSSA Method
    parameters_1 = ReversibleMichaelisMenten.Parameters(k_equilibrium=1.5)
    parameters_2 = ReversibleMichaelisMenten.Parameters(k_equilibrium=2.0)
    parameters_3 = ReversibleMichaelisMenten.Parameters(k_equilibrium=3.0)

    reaction1 = Reaction(name='E1',
                         mechanism=ReversibleMichaelisMenten,
                         reactants=metabolites_1,
                         )

    reaction2 = Reaction(name='E2',
                         mechanism=ReversibleMichaelisMenten,
                         reactants=metabolites_2,
                         )

    reaction3 = Reaction(name='E3',
                         mechanism=ReversibleMichaelisMenten,
                         reactants=metabolites_3,
                         )

    this_model = KineticModel()
    this_model.add_reaction(reaction1)
    this_model.add_reaction(reaction2)
    this_model.add_reaction(reaction3)

    the_boundary_condition = ConstantConcentration(this_model.reactants['A'])
    this_model.add_boundary_condition(the_boundary_condition)

    the_boundary_condition = ConstantConcentration(this_model.reactants['D'])
    this_model.add_boundary_condition(the_boundary_condition)

    this_model.parametrize_by_reaction({'E1': parameters_1,
                                        'E2': parameters_2,
                                        'E3': parameters_3})
    return this_model

@pytest.mark.dependency(name='build_linear_pathway_model')
def test_sampling_linear_pathway():
    this_model = build_linear_pathway_model()

    this_model.prepare(mca=True)
    this_model.compile_mca(sim_type = QSSA)


    flux_dict = {'E1': 1.0, 'E2': 1.0, 'E3': 1.0}
    concentration_dict = {'A': 3.0, 'B': 2.0, 'C': 1.0, 'D': 0.5}


    parameters = SimpleParameterSampler.Parameters(n_samples=10)
    sampler = SimpleParameterSampler(parameters)

    parameter_population_A = sampler.sample(this_model, flux_dict,
                                          concentration_dict, seed = 10)

    parameter_population_B = sampler.sample(this_model, flux_dict,
                                          concentration_dict, seed = 10)

    parameter_population_C = sampler.sample(this_model, flux_dict,
                                          concentration_dict, seed = 20)

    assert(parameter_population_A == parameter_population_B)
    assert( not(parameter_population_B == parameter_population_C))