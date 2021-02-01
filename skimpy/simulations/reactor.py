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
import yaml

from skimpy.io.yaml import load_yaml_model
from skimpy.core import BoundaryCondition
from skimpy.core.reactor import Reactor
from skimpy.utils.general import make_subclasses_dict
from skimpy.utils.tabdict import TabDict


def get_boundary_subclasses():
    return make_subclasses_dict(BoundaryCondition)


def make_batch_reactor(path):
    """

    :param path: path to config yml file to setup reactor
    :return:
    """
    with open(path,'r') as fid:
            the_dict = yaml.full_load(fid)

    # Load models
    models = dict()
    for name, file in the_dict['models'].items():
        model =  load_yaml_model(file)
        model.name = name
        models[name] = model

    # Load scaling
    concentration_scaling = float(the_dict['scaling']['concentration'])
    density = float(the_dict['scaling']['density'])
    gDW_gWW = float(the_dict['scaling']['gDW_gWW'])
    time_scaling = float(the_dict['scaling']['time'])
    flux_scaling =  1e-3 * (gDW_gWW * density) \
                    * concentration_scaling \
                    / time_scaling

    biomass_scaling = {k:flux_scaling/float(v) for k,v in
                       the_dict['biomass_scaling'].items()}

    biomass_reactions = dict()
    for name, biomass_rxn in the_dict['biomass'].items():
        biomass_reactions[name] = models[name].reactions[biomass_rxn]

    extracellular_compartment = the_dict['extracellular_compartment']
    reactor_volume = float(the_dict['reactor_volume'])

    reactor = Reactor(models.values(), biomass_reactions, biomass_scaling,
                      extracellular_compartment=extracellular_compartment)

    for model in reactor.models.values():
        model.compartments[extracellular_compartment].parameters.volume.value \
            = reactor_volume

    # Boundary Conditions
    for the_bc_dict in the_dict['boundary_conditions'].values():
        TheBoundaryCondition = get_boundary_subclasses()[the_bc_dict.pop('class')]
        reactant = reactor.variables[the_bc_dict.pop('reactant')]
        the_bc = TheBoundaryCondition(reactant, **the_bc_dict)
        reactor.add_boundary_condition(the_bc)

    # Init the reactor initial conditions in correct order
    reactor.initial_conditions = TabDict([(x,0.0) for x in reactor.variables])
    # Add medium:
    for met, conc in the_dict['initial_medium'].items():
        reactor.initial_conditions[met]  = float(conc) * concentration_scaling

    # Add scaling fator to the rectaor
    reactor.concentration_scaling = concentration_scaling
    reactor.flux_scaling = flux_scaling

    return reactor


