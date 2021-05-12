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

import ctypes
import re
import os

import numpy as np

import tempfile

import multiprocessing
from sympy.printing import ccode
from sympy import Symbol


"""
TODO: 
We currently use cython to generate the code, this offers a safe type conversion.
The generated cython code wraps every line in a function increasing the compilation
and execution time 
But as sympy code printer directly prints C-code we should use this to generate a swig interface.
"""

# This should be plat form depednent using distrtools
COMPILER = "gcc -fPIC -shared -w -O3"

# Test to write our own compiler
INCLUDE = "#include <stdlib.h>\n" \
          "#include <math.h>\n"

FUNCTION_DEFINITION_HEADER = "void function(double *input_array, double *output_array){ \n"
FUNCTION_DEFINITION_FOOTER = ";\n}"

def make_cython_function(symbols, expressions, quiet=True, simplify=True, optimize=False, pool=None):

    code_expressions = generate_vectorized_code(symbols,
                                                expressions,
                                                simplify=simplify,
                                                pool=pool)


    # Write the code to a temp file
    code = INCLUDE + FUNCTION_DEFINITION_HEADER + code_expressions + FUNCTION_DEFINITION_FOOTER
    path_to_c_file = write_code_to_tempfile(code)
    path_to_so_file = path_to_c_file.replace('.c', '.so')

    # Compile the code
    cmd = " ".join([COMPILER, '-o ',path_to_so_file, path_to_c_file] )
    # Todo catch errors
    # Todo catch errors
    os.system(cmd)

    # Import the function
    fun = ctypes.CDLL(path_to_so_file)

    def this_function(input_array,output_array):
        # Input pointers
        fun.function.argtypes = [ctypes.POINTER(ctypes.c_double),
                                 ctypes.POINTER(ctypes.c_double),]
        #Cast to numpy float
        if not type(input_array) ==  np.ndarray.dtype:
            input_array = np.array(input_array, dtype=np.float)

        #x.ctypes.data_as(ctypes.POINTER(ctypes.c_long))
        fun.function(input_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)) ,
                     output_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), )

    return this_function

def write_code_to_tempfile(code,file_path=None):
    if file_path is None:
        # make a tempfile
        (_, file_path) = tempfile.mkstemp(suffix = '.c')

    with open(file_path, "w") as text_file:
        text_file.write(code)
    return file_path

def generate_vectorized_code(inputs, expressions, simplify=True, optimize=False, pool=None):
    # input substitution dict:
    input_subs = {str(e): "input_array[{}]".format(i)
                  for i, e in enumerate(inputs)}

    if pool is None:
        cython_code = []
        for i,e in enumerate(expressions):
            if simplify:
                cython_code.append(generate_a_code_line_simplfied((i,e,input_subs)))
            else:
                cython_code.append(generate_a_code_line((i, e, input_subs)))

    else:
        if simplify:
            input_subs_input = [input_subs, ]*len(expressions)
            i,e =zip(*enumerate(expressions))
            cython_code = pool.map(generate_a_code_line_simplfied, zip(i,e,input_subs_input) )

        else:
            input_subs_input = [input_subs, ] * len(expressions)
            i, e = zip(*enumerate(expressions))
            cython_code = pool.map(generate_a_code_line, zip(i, e, input_subs_input))

    cython_code = ';\n'.join(cython_code)

    return cython_code


from sympy import cse

def generate_a_code_line_simplfied(input , optimize=False):
    i, e, input_subs = input

    # Use common sub expressions instead of simpilfy
    # Generate directly unique CSE Symbols and tranlate them to ccode
    if optimize:
        common_sub_expressions, main_expression = cse(e.simplify())
    else:
        common_sub_expressions, main_expression = cse(e)

    # Generate unique symbols for the common subexpressions
    cse_subs = {}
    common_sub_expressions_unique = []
    for this_cse in common_sub_expressions:
        gen_sym = str(this_cse[0])
        unique_sym = Symbol("cse_{}_{}".format(i,gen_sym))
        cse_subs[ this_cse[0]] = unique_sym
        common_sub_expressions_unique.append(
            [unique_sym,this_cse[1]] )

    # Substitute the cse symbols in mainexpression and other cse
    for this_cse in common_sub_expressions_unique:
        this_cse[1] = this_cse[1].subs(cse_subs)

    main_expression = main_expression[0].subs(cse_subs)

    cython_code = ''
    for this_cse in common_sub_expressions_unique:
        cython_code=cython_code+'double {} = {} ;\n'.format(str(this_cse[0]),
                                                    ccode(this_cse[1],standard='C99'))


    cython_code = cython_code+"output_array[{}] = {} ;".format(i,ccode(main_expression
                                                                      ,standard='C99')
                                                              )

    # Substitute integers in the cython code
    cython_code = re.sub(r"(\ |\+|[^e]\-|\*|\(|\)|\/|\,)([1-9])(\ |\+|\-|\*|\(|\)|\/|\,)",
                         r"\1 \2.0 \3 ",
                         cython_code)

    for str_sym, array_sym in input_subs.items():
        cython_code = re.sub(r"(\ |\+|\-|\*|\(|\)|\/|\,)({})(\ |\+|\-|\*|\(|\)|\/|\,)".format(str_sym),
                             r"\1 {} \3 ".format(array_sym),
                             cython_code)

    return cython_code


def generate_a_code_line(input, optimize=False):
    i, e, input_subs = input

    if optimize:
        cython_code = "output_array[{}] = {} ".format(i, ccode(e.simplify()), standard='C99')
    else:
        cython_code = "output_array[{}] = {} ".format(i,ccode(e, standard='C99'))


    # Substitute integers in the cython code
    cython_code = re.sub(r"(\ |\+|[^e]\-|\*|\(|\)|\/|\,)([1-9])(\ |\+|\-|\*|\(|\)|\/|\,)",
                         r"\1 \2.0 \3 ",
                         cython_code)

    for str_sym, array_sym in input_subs.items():
        cython_code = re.sub(r"(\ |\+|\-|\*|\(|\)|\/|\,)({})(\ |\+|\-|\*|\(|\)|\/|\,)".format(str_sym),
                             r"\1 {} \3 ".format(array_sym),
                             cython_code)
    return cython_code