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

import Cython
import re
import os

from sympy.printing import ccode

CYTHON_DECLARATION = "# cython: boundscheck=True, wraparound=False, language_level=3"+\
                     "nonecheck=True, initializecheck=False , optimize=False\n"+\
                     "# distutils: language = c++ \n"

SQRT_FUNCTION = "cdef extern from \"math.h\": \n double sqrt(double x) \n"

def _set_cflags():
    """ Suppress cython warnings by setting -w flag """

    flags = '-w -O0'

    if 'CFLAGS' not in os.environ:
        os.environ['CFLAGS'] = flags
    elif not flags in os.environ['CFLAGS']:
        os.environ['CFLAGS'] = os.environ['CFLAGS'] + " " + flags



def make_cython_function(symbols, expressions, quiet=True , simplify=True):

    code_expressions = generate_vectorized_code(symbols, expressions, simplify=simplify)

    _set_cflags()

    def this_function(input_array,output_array):
        Cython.inline(CYTHON_DECLARATION\
                      +SQRT_FUNCTION\
                      +code_expressions,
                      quiet=quiet)

    return this_function


def generate_vectorized_code(inputs, expressions, simplify=True):
    # input substitution dict:
    input_subs = {str(e):"input_array[{}]".format(i)
                  for i,e in enumerate(inputs)}

    if simplify:
        cython_code = ["output_array[{}] = {} ".format(i,ccode(e.simplify()))
                       for i,e in enumerate(expressions)]
    else:
        cython_code = ["output_array[{}] = {} ".format(i, ccode(e))
                       for i, e in enumerate(expressions)]

    cython_code = '\n'.join(cython_code)

    for str_sym, array_sym in input_subs.items():
        cython_code = re.sub(r"(\ |\+|\-|\*|\(|\)|\/|\,)({})(\ |\+|\-|\*|\(|\)|\/|\,)".format(str_sym),
                             r"\1 {} \3 ".format(array_sym),
                             cython_code)
    # Substitute integers in the cython code
    cython_code = re.sub(r"(\ |\+|\-|\*|\(|\)|\/|\,)([1-9])(\ |\+|\-|\*|\(|\)|\/|\,)",
                         r"\1 \2.0 \3 ",
                         cython_code)
    return cython_code
