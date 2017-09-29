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
def uni_uni_from_thermo(name,substrates,thermo_data):
    """Create an reaction from thermodynamic data"""

    # v_maxes
    v_max_f =  thermo_data['flux'] \
                /((1.0-thermo_data['gamma'])*thermo_data['sig_S'])

    v_max_r =  thermo_data['gamma']*thermo_data['flux'] \
                /((1.0-thermo_data['gamma'])*thermo_data['sig_P'])
    # K_M values
    nominator = (1.0-thermo_data['sig_S']-thermo_data['sig_P'])
    K_S = thermo_data['S']*nominator/thermo_data['sig_S']
    K_P = thermo_data['P']*nominator/thermo_data['sig_P']

    params = {'K_S'     : K_S ,
              'K_P'     : K_P ,
              'v_max_r' :v_max_r,
              'v_max_f' :v_max_f,
              'E_tot'   :thermo_data['E_tot']}

    new_enzyme = ReversibleMichaelisMenten(name,substrates,params)
    return new_enzyme

def uni_uni_from_data(name,substrates,params):
    """Create an enzmye from enzyme data"""
    new_enzyme = ReversibleMichaelisMenten(name,substrates,params)

    return new_enzyme