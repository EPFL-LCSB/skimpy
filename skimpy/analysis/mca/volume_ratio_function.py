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

from numpy import array, double, zeros
from numpy import append as append_array

from sympy import symbols

from skimpy.utils.tabdict import TabDict
from skimpy.utils.compile_sympy import make_cython_function

class VolumeRatioFunction:
    def __init__(self, model, variables,  parameters, pool=None):
        """
        Constructor for a precompiled function to compute elasticities
        numerically
        :param variables: a list of strings denoting
                                      the independent variables names
        :param expressions: dict of  non-zero sympy expressions for the rate of
                            change of a variable indexed by a tuple of the matrix position
                            e.g: (1,1)
        :param parameters:  list of parameter names
        :param shape: Tuple defining the over all matrix size e.g (10,30)

        """
        self.reactants = model.reactants

    def __call__(self,parameters):
        """
        Return a list of volume ratios
        """
        values = [parameters[v.compartment.parameters.cell_volume.symbol] /
                  parameters[v.compartment.parameters.volume.symbol]
                  for k, v in self.reactants.items()]

        return array(values)

