from skimpy.nullspace import left_integer_nullspace

import numpy as np
from scipy.sparse import random
from scipy import stats

class CustomRandomState(object):
    def randint(self, k):
        i = np.random.randint(k)
        return i

rs = CustomRandomState()
rvs = stats.poisson(2, loc=10).rvs
S = random(5,6, density=0.1, random_state=rs, data_rvs=rvs)

print(S.todense())

ns = left_integer_nullspace(S.todense())

print(ns)
null = ns @ S.todense()

print(null)
print(np.any(null))