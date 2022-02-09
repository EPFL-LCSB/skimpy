"""


"""

from skimpy.core.parameters import ParameterValuePopulation
from scipy.stats import multivariate_normal

import tensorflow as tf

import pandas as pd
import numpy as np

EPSILON = 1e-9

class SecureMultivariateNormal(object):
    def __init__(self, mu, sigma, var):
        self.variable_parameters = var > EPSILON
        self.constant_parameters = var < EPSILON

        self.mu = mu

        self.variable_index = var.index[self.variable_parameters]
        self.const_index = var.index[self.constant_parameters]

        # TODO: RAISE A MORE INFORMATIVE ERROR WHEN THE COV IS SIGULAR!

        self._dist = multivariate_normal(mu[self.variable_parameters],
                                         sigma.loc[self.variable_parameters,self.variable_parameters])

    def rvs(self,size,random_state=None):
        values_var = self._dist.rvs(size=size, random_state=random_state)
        if size > 1:
            df = pd.DataFrame(values_var, columns=self.variable_index)
        else:
            df = pd.DataFrame(values_var, index=self.variable_index, columns=[0]).T
        values_cons = pd.concat( [ self.mu[self.constant_parameters] , ]*size , axis=1)

        return pd.concat([df, values_cons.T], axis=1)


class LogNormalPriorParameterDistribution():
    """
    TF Based model
    """
    pass

class PosteriorLogNormalParameterPopulation(object):
    """

    """
    def __init__(self, parameter_poulations, likelyhoods=None):
        self.mu = []
        self.sigma = []
        self.pdf = []

        for pop in parameter_poulations:
            var = pop.log_var()
            mu = pop.log_mean()
            sigma = pop.log_cov()

            self.mu.append( mu )
            self.sigma.append( sigma )
            self.pdf.append( SecureMultivariateNormal(mu,sigma,var) )

        N = len(parameter_poulations)
        if likelyhoods is None:
            self.weights = np.ones(N)/N
        else:
            self.weights = likelyhoods/sum(likelyhoods)

        self.cum_weights = np.cumsum(self.weights)

    def resample(self, N, seed=None):
        if not seed is None:
            np.random.seed(seed=seed)

        # Gillespie type of algorithm for resampling
        x = []
        for s in range(N):
            r = np.random.rand()
            j = sum(r > self.cum_weights)
            x_i = np.exp( self.pdf[j].rvs(size=1, random_state=None) )
            x.append(x_i)
        df = pd.concat(x, axis=0, ignore_index=True, )
        return df



class PosteriorNormalParameterPopulation(object):
    """

    """
    def __init__(self, parameter_poulations, likelyhoods=None):
        self.mu = []
        self.sigma = []
        self.pdf = []

        for pop in parameter_poulations:
            var = pop.var()
            mu = pop.mean()
            sigma = pop.cov()

            self.mu.append( mu )
            self.sigma.append( sigma )
            self.pdf.append( SecureMultivariateNormal(mu,sigma,var) )

        N = len(parameter_poulations)
        if likelyhoods is None:
            self.weights = np.ones(N)/N
        else:
            self.weights = likelyhoods/sum(likelyhoods)

        self.cum_weights = np.cumsum(self.weights)

    def resample(self, N, seed=None):
        if not seed is None:
            np.random.seed(seed=seed)

        # Gillespie type of algorithm for resampling
        x = []
        for s in range(N):
            r = np.random.rand()
            j = sum(r > self.cum_weights)
            x_i = self.pdf[j].rvs(size=1, random_state=None)
            x.append(x_i)
        df = pd.concat(x, axis=0, ignore_index=True, )
        return df