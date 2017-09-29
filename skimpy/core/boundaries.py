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

class Boundary(object):
    def __init__(self,boundary_method):
        self.boundary_method = boundary_method

    def __call__(self,expressions):
        self.boundary_method(expressions)


class ConstantConcentration(Boundary):

    def __init__(self,metabolite):
        method = lambda expression: set_const(expression,metabolite)
        Boundary.__init__(self,method)


# Private function
def set_constant_concentration(expressions,metabolite):
    expressions[metabolite] = expressions[metabolite]*0.0
