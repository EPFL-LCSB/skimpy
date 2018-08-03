from skimpy.nullspace import integer_nullspace

import numpy as np


A = np.random.randint(0,10,size=(3,3))

A = np.array([[1,0,0],[0,0,1],[1,0,1]])

print(A)

ns = integer_nullspace(A)

print(ns)

print(A.dot(ns))