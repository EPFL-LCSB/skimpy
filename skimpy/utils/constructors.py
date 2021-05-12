# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

Reaction constructors, utility funs

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

from ..mechanisms import *

#Create enzyme parameterization from
#Thermodynamics data
def reversible_michaelis_menten_from_thermo_data(thermo_data):
    """Create reaction parameters from thermodynamic data"""

    # v_maxes
    vmax_forward =  thermo_data['flux'] \
                    /((1.0-thermo_data['thermo_displacement'])\
                    *thermo_data['saturation_substrate'])

    vmax_backward =  thermo_data['thermo_displacement']*thermo_data['flux'] \
                    /((1.0-thermo_data['thermo_displacement'])\
                    *thermo_data['saturation_substrate'])
    # K_M values
    nominator = 1.0-thermo_data['saturation_product'] \
                        -thermo_data['saturation_substrate']

    km_substrate = thermo_data['concetration_substrate']\
                        *nominator/thermo_data['saturation_substrate']

    km_product = thermo_data['concetration_product']\
                        *nominator/thermo_data['saturation_product']


    parameters = ReversibleMichaelisMenten.Parameters(
        vmax_forward = vmax_forward,
        vmax_backward = vmax_backward,
        km_substrate = km_substrate,
        km_product = km_product,
        total_enzyme_concentration = thermo_data['total_enzyme_concentration'],
    )


    return parameters
