# -*- coding: utf-8 -*-
"""
.. module:: pytfa
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


import numpy as np

from scipy.sparse import eye as speye
from scipy.sparse import find as sfind
from scipy.sparse import hstack, vstack
from skimpy.nullspace import left_integer_nullspace
from sympy import Matrix, nsimplify

def rational_left_basis(S):
    S = np.array(S, dtype=np.int)
    return integer_nullspace(S.T)


"""
Disclaimer:

This is an ugly recoding of MATLAB's sbioconsmoiety.m, 'semipos' method
"""

def integer_left_basis_semipos(S):

    #TODO: chack sparsity

    n,r = S.shape

    T = hstack((speye(n), S),format='lil')

    for j in range(n+r-1,n-1,-1):

        zero_rows,_ = T[:,j].nonzero()
        zero_ix = list(set(range(T.shape[0])).difference(zero_rows))
        Tnew = T[zero_ix, 0:j]
        posinds,_,_ = sfind(T[:,j]>0)
        neginds,_,_ = sfind(T[:,j]<0)

        lni = len(neginds)

        for i in posinds:
            nnz_rows, nnz_cols = (T[[i]*lni,:n+1] + T[neginds,:n+1]).nonzero()
            zero_cols = list(set(range(T.shape[1])).difference(nnz_cols))

            for k in range(lni):
                flags = T[:, zero_cols[k]]
                flags[i] = True
                flags[neginds[k]] = True

                if flags.nnz == flags.shape[0]:
                    newrow = -1*T[neginds[k],j]*T[i,:j] + T[i,j]*T[neginds[k],:j]
                    Tnew = vstack(Tnew, newrow)

        T = Tnew

    g=np.zeros((T.shape[0],1))

    for c in range(len(g)):
        g[c] = local_gcd(T[c,:])

    return T/g,T,g



def local_gcd(v):
    # if not v == np.round(v):
    #     return 1

    x = v[v.nonzero()]
    g=0
    for xi in x.transpose():
        g = numpy_gcd(g,xi.tolist()[0])
        if g==1:
            return g

    return g

def numpy_gcd(a, b):
    """
    https://stackoverflow.com/questions/15569429/numpy-gcd-function

    :param a:
    :param b:
    :return:
    """
    a_, b_ = np.broadcast_arrays(a, b)
    a_ = a_.copy()
    b_ = b_.copy()
    pos = np.nonzero(b_)[0]
    while len(pos) > 0:
        b2 = b_[pos]
        a_[pos], b_[pos] = b2, a_[pos] % b2
        pos = pos[b_[pos]!=0]
    return a_