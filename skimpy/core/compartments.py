# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2020 Laboratory of Computational Systems Biotechnology (LCSB),
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
from skimpy.core.itemsets import Parameter
from skimpy.utils import TabDict, iterable_to_tabdict

class Compartment(object):
    """
    Global Compartment class
    """
    def __init__(self, name, model=None ):

        #TODO Do we need additional info for the future?
        self.name = name

        #Parameters so far only volume but maybe others soon
        volume = Parameter(name='volume', model=model, value=None, suffix=name  )
        cell_volume = Parameter(name='cell_volume', model=model, value=None, suffix=name )

        self.parameters = iterable_to_tabdict([volume, cell_volume])

