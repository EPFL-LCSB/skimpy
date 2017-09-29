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

from scipy.integrate import ode
from numpy import array, append

from .ode_fun import make_ode_fun
from .solution import Solution


class KineticModel(object):
    # Consult with Pierre about this class!
    # Better use dicts with names! + inherited objects!

    def __init__(self, metabolites = [], reactions = [], boundaries = []):
        # initialize the model is stated
        # FIXME Add dictlists from cobra ? or reimplement a similar data structure
        self.metabolites = metabolites    #List of metabolite objects/ids
        self.reactions   = reactions      #List of enzyme objects
        self.boundaries  = boundaries
        self._modifed    = True

    def add_reaction(self, reaction):
        # Add an enzyme to the model
        if has_same_name(reaction,self.reactions):
            reaction.name = reaction.name+"_double"

        self.reactions.append(reaction)
        for this_metabolite in reaction.metabolites:
            self.metabolites.append(this_metabolite)
        self._modifed = True

    def add_boundary(self, boundary):
        # Add a boundary to the model
        self.boundaries.append(boundary)

    def solve_ode(self,
                  time_int,
                  initial_concentrations,
                  sim_type = 'QSSA',
                  solver_type = 'vode',
                  reltol = 1e-8,
                  abstol = 1e-8):

        # Create the ode_fun if modified or non exisitent
        if self._modifed and not(hasattr(self,'ode_fun')):
            self.ode_fun = make_ode_fun(self,sim_type)
            self._modifed = False

        # Choose a solver
        self.solver = get_ode_solver(self.ode_fun,solver_type,reltol,abstol)

        # solve the ode
        t_sol,y_sol = solve_ode(self.solver,time_int,initial_concentrations)

        return Solution(self,t_sol,y_sol)


def get_ode_solver(  ode_fun,
                     solver_type = "vode",
                     reltol = 1e-8,
                     abstol = 1e-8):

    # Initilize the integrator
    ode_solver = ode(ode_fun)
    # Set proprties
    ode_solver.set_integrator(  solver_type,
                                method='bdf',
                                atol=abstol,
                                rtol=reltol )

    return ode_solver



def solve_ode(solver,time_int,initial_concentrations):

    solver.set_initial_value(initial_concentrations, time_int[0])

    t_sol = [time_int[0]]
    y_sol = [initial_concentrations]

    while solver.t <= time_int[1] and solver.successful():
            solver.integrate(time_int[1], step=True)
            t_sol.append(solver.t)
            y_sol.append(solver.y)

    return t_sol,y_sol



def has_same_name(reaction,reactions):
    for this_reaction in reactions:
        if this_reaction.name == reaction.name:
            return True
    return False