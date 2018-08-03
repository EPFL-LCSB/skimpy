
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

    cdef slong fmpz_mat_nullspace(fmpz_mat_t N, fmpz_mat_t A);

    cdef void fmpz_mat_init(fmpz_mat_t A, long rows, long columns);
    cdef fmpz_t fmpz_mat_entry(fmpz_mat_t A, long row, long column);
    cdef void fmpz_set_ui(fmpz_t A, slong x);
    cdef slong fmpz_get_ui(fmpz_t A);
    cdef void fmpz_mat_clear(fmpz_mat_t A);


cdef int get_fmpz_mat_entry(fmpz_mat_t A, long row, long column):
    cdef fmpz_t entry_t = fmpz_mat_entry( A, row, column);
    cdef int entry = fmpz_get_ui(entry_t );
    return entry


cpdef integer_nullspace(matrix):

    n = matrix.shape[0]
    m = matrix.shape[1]

    cdef fmpz_mat_t matrix_t;
    cdef fmpz_mat_t nullspace_t;

    fmpz_mat_init(matrix_t, n, m)
    fmpz_mat_init(nullspace_t, n, m)

    for i in range(n):
        for j in range(m):
            entry = fmpz_mat_entry(matrix_t, i, j)
            fmpz_set_ui(entry , matrix[i,j])


    cdef rank_def = fmpz_mat_nullspace(nullspace_t, matrix_t)


    nullspace = np.zeros( (n,rank_def) )
    for i in range(n):
        for j in range(rank_def):
            nullspace[i,j] = get_fmpz_mat_entry(nullspace_t,i,j)

    #Clear matrices
    fmpz_mat_clear(nullspace_t)
    fmpz_mat_clear(matrix_t)

    return nullspace


