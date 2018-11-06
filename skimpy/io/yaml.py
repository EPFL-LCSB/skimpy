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

import yaml


FIELDS_TO_SERIALIZE = ['variables',
                       # 'ode_fun',
                       'boundary_conditions',
                       'parameters',
                       #'logger',
                       # '_simtype',
                       'constraints',
                       'solver',
                       # '_recompiled', '_modifed',
                       'reactions',
                       # '_modified',
                       'name', 'initial_conditions']


def export_to_yaml(model, path=None, **kwargs):

    dict_model = vars(model)
    fields_not_to_serialize = [x for x in dict_model if not x in FIELDS_TO_SERIALIZE]
    [dict_model.pop(k) for k in fields_not_to_serialize]

    if path is not None:
        with open(path, 'w') as fid:
            yaml.dump(dict_model, fid)
            return True
    else:
        return yaml.dump(dict_model)

def load_yaml_model(path):
    pass