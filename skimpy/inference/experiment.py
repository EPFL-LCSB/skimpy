"""

"""

from skimpy.utils.tabdict import TabDict, iterable_to_tabdict
from scipy.stats import multivariate_normal

import numpy as np
import pandas as pd

class Observable(object):
    def __init__(self, name, dict):
        self.name = name
        self._dict = dict

    def __call__(self, X):
        return sum([X[k]*v for k,v in self._dict.items() ])


class SteadyStateExperiment(object):
    def __init__(self,observables,mu,var):

        self.observables = iterable_to_tabdict(observables)
        self.mu = mu
        self.var = var

        self._distribution_observables = multivariate_normal(self.mu,
                                                             np.diag(self.var))

    def p(self,fluxes,concentrations):

        #TODO: CORRECT FOR DIFFERENT SAMPLE SIZES!

        X_i = pd.concat([fluxes, concentrations], axis=1)
        n_obs = len( self.observables)
        n_samples = X_i.shape[0]
        observables_i = np.zeros( (n_obs, n_samples))

        for i, observable in enumerate(self.observables.values()):
            observables_i[i,:] = observable(X_i)


        p_i = self._distribution_observables.pdf(observables_i.T)

        return p_i

