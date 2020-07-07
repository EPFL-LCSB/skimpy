
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
import numpy as np

cdef extern from "nullspace.pyh":
    ctypedef long slong
    ctypedef slong fmpz
    ctypedef fmpz fmpz_t[1]

    ctypedef struct fmpz_mat_struct:
        fmpz * entries;
        slong r;
        slong c;
        fmpz ** rows;

    ctypedef fmpz_mat_struct fmpz_mat_t[1]

    # Rationals
    ctypedef struct fmpq:
         fmpz num;
         fmpz den;

    ctypedef fmpq fmpq_t[1]

    ctypedef struct fmpq_mat_struct:
        fmpq * entries;
        slong r;
        slong c;
        fmpq ** rows;

    ctypedef fmpq_mat_struct fmpq_mat_t[1]

    cdef void fmpz_mat_init(fmpz_mat_t A, int rows, int columns);
    cdef fmpz_t fmpz_mat_entry(fmpz_mat_t A, int row, int column);
    cdef void fmpz_set_ui(fmpz_t A, slong x);
    cdef slong fmpz_get_ui(fmpz_t A);
    cdef void fmpz_mat_clear(fmpz_mat_t A);

    cdef void fmpz_clear(fmpz_t x)
    cdef void fmpz_init(fmpz_t x)
    cdef void fmpz_set(fmpz_t dest , const fmpz_t src)


    cdef void fmpq_mat_init(fmpq_mat_t A, int rows, int columns);
    cdef fmpq_t fmpq_mat_entry(fmpq_mat_t A, int row, int column);
    cdef void fmpq_init(fmpq_t x)
    cdef void fmpq_set(fmpq_t dest , const fmpq_t src)
    cdef void fmpq_clear(fmpq_t x)
    cdef void fmpq_mat_clear(fmpq_mat_t A);

    cdef slong fmpq_mat_rref(fmpq_mat_t B, const fmpq_mat_t A)

    cdef slong fmpz_mat_nullspace ( fmpz_mat_t B , const fmpz_mat_t A )



cdef int get_fmpz_mat_entry(fmpz_mat_t A, int row, int column):
    cdef fmpz_t entry_t = fmpz_mat_entry( A, row, column);
    cdef int entry = fmpz_get_ui(entry_t );
    return entry


cpdef left_integer_nullspace(matrix):

    if not np.issubdtype(matrix.dtype,np.integer):
        raise TypeError("The Matrix is not integer.")

    n_0 = matrix.shape[0]
    m_0 = matrix.shape[1]

    matrix = np.concatenate((matrix,np.eye(n_0)) ,axis=1)

    n = matrix.shape[0]
    m = matrix.shape[1]

    cdef fmpq_mat_t matrix_t;
    cdef fmpq_mat_t rrechelon_t;


    fmpq_mat_init(matrix_t, n, m)
    t = max([m,n])
    fmpq_mat_init(rrechelon_t, n, m)

    cdef fmpq_t q
    fmpq_init(q)

    for i in range(n):
        for j in range(m):
            entry = fmpq_mat_entry(matrix_t, i, j)
            q.num = matrix[i,j]
            q.den = 1
            fmpq_set(entry , q )

    cdef slong* perm
    cdef fmpz_t den
    cdef rank = fmpq_mat_rref(rrechelon_t, matrix_t)

    echelon = np.zeros( (n,m) )

    for i in range(n):
        for j in range(m):
            q = fmpq_mat_entry(rrechelon_t,i,j)
            echelon[i,j] = float(q.num)/float(q.den)


    fmpq_mat_clear(rrechelon_t)
    fmpq_mat_clear(matrix_t)
    fmpq_clear(q)

    rref = echelon[:,:m_0]
    transformation = echelon[:,m_0:]

    zero_rows = [i for i,row in enumerate(rref) if not np.any(row)]

    left_nullspace = transformation[zero_rows,:]

    return left_nullspace


cpdef right_integer_nullspace(matrix):

    if not np.issubdtype(matrix.dtype,np.integer):
        raise TypeError("The Matrix is not integer.")

    n = matrix.shape[0]
    m = matrix.shape[1]


    cdef fmpz_mat_t matrix_t;
    cdef fmpz_mat_t nullspace_t;


    fmpz_mat_init(matrix_t, n, m)
    fmpz_mat_init(nullspace_t, m, m)

    cdef fmpz_t z
    fmpz_init(z)

    for i in range(n):
        for j in range(m):
            fmpz_set_ui ( fmpz_mat_entry (matrix_t , i , j ) , matrix[i,j] )

    cdef rank = fmpz_mat_nullspace(nullspace_t, matrix_t)

    nullspace = np.zeros( (rank,m) )

    for i in range(rank):
        for j in range(m):
            z = fmpz_mat_entry(nullspace_t,i,j)
            value =fmpz_get_ui(z)
            nullspace[i,j] = value

    fmpz_mat_clear(nullspace_t)
    fmpz_mat_clear(matrix_t)
    fmpz_clear(z)

    return nullspace