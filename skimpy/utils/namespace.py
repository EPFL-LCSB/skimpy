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
""" Simulation types """
QSSA = 'qssa'
TQSSA = 'tqssa'
MCA  = 'mca'
ODE  = 'ode'
ELEMENTARY = 'elementary'


""" Jacobian Types"""
NUMERICAL = 'numerical'
SYMBOLIC = 'symbolic'

""" MCA Types """
NET = 'net'
SPLIT = 'split'


""" Item types """
PARAMETER = 'parameter'
VARIABLE  = 'variable'

""" Units """
KCAL = 'kcal'
KJ   = 'kJ'
JOULE = 'JOULE'


""" OTHER """
WATER_FORMULA = 'H2O'

