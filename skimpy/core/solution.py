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
from ..viz.plotting import timetrace_plot, plot_population_per_variable
from ..utils import TabDict,iterable_to_tabdict

from copy import deepcopy

import pandas as pd

# Class for ode solutions
class ODESolution:
    def __init__(self, model, solution):
        self.ode_solution = solution

        self.time    = np.array(solution.values.t)

        self.species = np.array(solution.values.y)
        self.names = [x for x in model.ode_fun.variables]

        # TODO: Cleanup this
        concentrations = iterable_to_tabdict([])
        for this_species, this_name in zip(self.species.T, self.names):
            concentrations[this_name] = this_species

        self.concentrations = pd.DataFrame.from_dict(concentrations, orient='columns')

    def plot(self, filename='', **kwargs):
        timetrace_plot(self.time, self.species, filename, legend=self.names, **kwargs)

    def copy(self):
        return deepcopy(self)


class ODESolutionPopulation:

    def __init__(self, list_of_solutions, index=None):

        sol_cols = list(list_of_solutions[0].concentrations.keys())
        self.data = pd.DataFrame(columns=['solution_id', 'time']+sol_cols)

        temp = []
        if index is None:
            for e, td in enumerate(list_of_solutions):
                new_block = pd.DataFrame.from_dict(td.concentrations.copy())
                new_block['time'] = td.time
                new_block['solution_id'] = e
                temp.append(new_block)

        else:
            for e, td in zip(index,list_of_solutions):
                new_block = pd.DataFrame.from_dict(td.concentrations.copy())
                new_block['time'] = td.time
                new_block['solution_id'] = e
                temp.append(new_block)

        self.data = pd.concat(temp, axis = 0)

    def plot(self, filename, variables=None, **kwargs):
        plot_population_per_variable(self.data, filename, variables=variables, **kwargs)

