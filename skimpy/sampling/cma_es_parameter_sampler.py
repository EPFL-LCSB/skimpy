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
WITHOUT WARRANTIE CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""
from collections import namedtuple
import numpy as np
#from scipy.sparse.linalg import eigs as eigenvalues
from scipy.linalg import eigvals as eigenvalues
from sympy import sympify, Symbol

from skimpy.sampling.utils import calc_max_eigenvalue, calc_parameters
from skimpy.utils.namespace import *

import random, array
from deap import algorithms



from skimpy.sampling import ParameterSampler, SaturationParameterFunction, FluxParameterFunction


class CMAESParameterSampler(ParameterSampler):
    """
    A simple parameter sampler that samples stable model parameters
    with respect to a steady state flux and concentration state
    """

    Parameters = namedtuple('Parameters', ['n_samples'])
    # TODO Talk to Pierre / Misko about simple sampler parameters
    # if parameters are not defined put default values
    Parameters.__new__.__defaults__ = (None,) * len(Parameters._fields)

    def sample(self,
               compiled_model,
               flux_dict,
               concentration_dict,
               seed=123,
               max_generation=10,
               sigma = 0.1,
               lambda_ = 1000,
               nhof = 100,
               max_eigenvalue = 0,
               min_km = 1e-3,
               max_km = 1e3,
               ):

        """

        :param compiled_model:
        :param flux_dict:
        :param concentration_dict:
        :param seed:
        :param max_generation:
        :param mutation_probability:
        :param eta:
        :return:
        """
        #
        from deap import base
        from deap import creator
        from deap import tools
        from deap import cma

        self.seed = seed
        random.seed(self.seed)

        symbolic_concentrations_dict = {Symbol(k):v
                                        for k,v in concentration_dict.items()}

        #Compile functions
        self._compile_sampling_functions(
            compiled_model,
            symbolic_concentrations_dict,
            flux_dict)

        """
        """

        self.compiled_model = compiled_model
        self.concentration_dict = concentration_dict
        self.flux_dict= flux_dict

        self.max_eigenvalue = max_eigenvalue

        """
        Define the optimization problem directly on the parameters 
        """

        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
        creator.create("Individual", array.array, typecode='d', fitness=creator.FitnessMin)

        self.sampling_parameters = [v for k,v in compiled_model.parameters.items()
                                     if (v.hook is not None) and (v.value is None) ]

        n_dim = len(self.sampling_parameters)

        bound_low = [min_km,]*n_dim
        bound_up  = [max_km,]*n_dim

        # Transform the bounds the logspace
        for i, the_parameter in enumerate(self.sampling_parameters):
            lb,ub = compiled_model.parameters[str(the_parameter.symbol)].bounds

            bound_low[i] = np.log(concentration_dict[the_parameter.hook.name]/ub) \
                if ub is not None else np.log(bound_low[i])

            bound_up[i] = np.log(concentration_dict[the_parameter.hook.name]/lb) \
                if lb is not None else np.log(bound_up[i])

        toolbox = base.Toolbox()
        toolbox.register("attr_float", init_parameters, bound_low, bound_up)
        toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.attr_float)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        #if hasattr(compiled_model,'pool'):
        #    toolbox.register("map", compiled_model.pool.map)

        toolbox.register("evaluate", self.fitness)

        parent = toolbox.individual()
        toolbox.evaluate(parent)

        strategy = cma.StrategyOnePlusLambda(parent=parent, sigma=sigma, lambda_=lambda_ )

        toolbox.register("generate", strategy.generate, creator.Individual)
        toolbox.register("update", strategy.update)

        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", np.nanmean)
        stats.register("std", np.nanstd)
        stats.register("min", np.nanmin)
        stats.register("max", np.nanmax)

        hof = tools.HallOfFame(nhof)

        result_parameters, _ = run_ea(toolbox, max_generation, stats=stats, hof=hof, verbose=True)

        parameter_population = []

        #TODO prune parameters sets that dont give eigenvalues

        for this_parameters in result_parameters:
            parameter_population.append(self.update_parameters(this_parameters))
        return parameter_population

    # Under construction new sampling with compiled function
    def _compile_sampling_functions(self,model,
                                    concentrations,
                                    fluxes):
        """
        Compliles the function for sampling using theano
        :param model:
        """

        model.flux_parameter_function = FluxParameterFunction(model,
                                                              model.parameters,
                                                              concentrations,)

    def fitness(self, parameters):

        # Get all parameters
        parameter_sample = self.update_parameters(parameters)

        lambda_max = calc_max_eigenvalue(parameter_sample,
                                         self.compiled_model,
                                         self.concentration_dict,
                                         self.flux_dict)

        if lambda_max < self.max_eigenvalue :
            return (self.max_eigenvalue,)
        else :
            return (lambda_max,)


    def update_parameters(self, parameters):
        parameter_sample = {v.symbol: v.value for k, v in self.compiled_model.parameters.items()}

        # Set Km parameter values (fitting in relative log1(Km/S) log space)
        for p, v in zip(self.sampling_parameters, parameters):
            parameter_sample[p.symbol] = self.concentration_dict[p.hook.name] / np.exp(v)

        symbolic_concentrations_dict = {Symbol(k): v
                                        for k, v in self.concentration_dict.items()}
        # Update the concentrations which are parameters (Boundaries)
        for k, v in symbolic_concentrations_dict.items():
            parameter_sample[k] = v

        for this_reaction in self.compiled_model.reactions.values():
            vmax_param = this_reaction.parameters.vmax_forward
            parameter_sample[vmax_param.symbol] = 1

        # Calculate the Vmax's
        self.compiled_model.flux_parameter_function(
            self.compiled_model,
            parameter_sample,
            symbolic_concentrations_dict,
            self.flux_dict
        )
        return parameter_sample


"""
Utils
"""


def run_ea(toolbox, ngen=None ,stats=None, hof=None, verbose=False):
    return algorithms.eaGenerateUpdate(toolbox, ngen=ngen, stats=stats, halloffame=hof)



"""
From DEAP tutorial 
"""
def init_parameters(low, up):
       return [random.uniform(a, b) for a, b in zip(low, up)]

def pareto_dominance(x,y):
    return tools.emo.isDominated(x.fitness.values, y.fitness.values)
