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
from .namespace import *
from .constants import *

def deltag0_to_keq(deltag0,temp,unit=KCAL,gas_constant=None):
    if gas_constant is not None:
        return np.exp(-1*deltag0/temp/gas_constant)
    if unit is KCAL:
        return np.exp(-1*deltag0*KCAL_IN_JOULE\
                      /GENERAL_GAS_CONSTANT/temp)
    if unit is KJ:
        return np.exp(-1*deltag0 * KJ_IN_JOULE\
                      /GENERAL_GAS_CONSTANT / temp)
    if unit is J:
        return np.exp(-1*deltag0 \
                      /GENERAL_GAS_CONSTANT / temp)
    else:
        raise ValueError
