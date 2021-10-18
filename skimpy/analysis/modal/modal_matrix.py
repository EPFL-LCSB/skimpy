# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2018 Laboratory of Computational Systems Biotechnology (LCSB),
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

from pandas import DataFrame, Index
import numpy as np
import scipy

from numpy.linalg import eig
from skimpy.utils import TabDict
from skimpy.analysis.ode import FluxFunction

def modal_matrix(kmodel,concentration_dict,parameters, flux_modes=False, tolerance=1e-12):
    """
    This function computes the transformation matrix W as described in Chapter 4 of
    Heinrich, Reinhart, and Stefan Schuster. The regulation of cellular systems.
    Springer Science & Business Media, 2012.

    The matrix W describes the composition of each pool variable X[i] whos
    dynamics are given by the eigenvalue lam[i] as a linear combination of
    the model variables S, which are concentrations.


    :param kmodel: A skimpy.core.KineticModel with compiled mca functions
    :param concentration_dict: A dict {'variable_name': value,  }
    :param parameters: A dict {'param_name': value,  }
    :return:
    """
    if not hasattr(kmodel,'jacobian_fun'):
        raise RuntimeError("MCA function not compiled cannot proceed modal analysis!")

    if not hasattr(kmodel,'flux_fun'):

        flux_expressions = TabDict([ (rxn.name, rxn.mechanism.reaction_rates['v_net'])
                             for rxn in kmodel.reactions.values()])
        flux_parameters = TabDict([(p.name, p) for rxn in kmodel.reactions.values()
                            for p in rxn.mechanism.expression_parameters])
        kmodel.flux_fun = FluxFunction(kmodel.variables,
                                       flux_expressions,
                                       flux_parameters,
                                       kmodel.pool)



    fluxes = kmodel.flux_fun(concentration_dict,parameters=parameters)

    concentrations= [concentration_dict[str(k)] for k in kmodel.variables]

    # Sort flux according to reactions
    fluxes = [fluxes[rxn] for rxn in kmodel.reactions]
    jacobian = kmodel.jacobian_fun(fluxes, concentrations, parameters, flux_jacobian=flux_modes)

    # The transformation matrix W is given by the eigenrows of the jacobian M,
    # or by the eigenvalues of the transpose of the jacobian M.

    # lam_1,_ = eig(jacobian.todense())
    # _  ,W_1 = eig(jacobian.todense().T)

    # # Check for small rows:
    # MAX_JAC = np.max(np.abs(jacobian.todense()))
    #
    # non_zero_rows = np.where(~np.all((np.abs(jacobian.todense()) < MAX_JAC*tolerance), axis=1))[0]
    #
    # reduced_jacobian = jacobian[non_zero_rows, :][:,non_zero_rows]

    lamr, lami, vl, vr, info = scipy.linalg.lapack.dgeev(jacobian.todense(),
                                                         compute_vl=1,
                                                         compute_vr=0,
                                                         )

    # Calculate the complex eigenvalues
    lam = lamr + 1j * lami

    # Calculate the comple eigenvetors
    W = np.zeros(vl.shape, dtype=complex)

    i = 0

    while i < len(lami):
        if lami[i] !=0:
            W[:, i]   = vl[:, i] + vl[:, i+1]*1j
            W[:, i+1] = vl[:, i] - vl[:, i+1]*1j
            i += 2
        else:
            W[:, i] = vl[:, i] + 0j
            i += 1

    W = W.T

    index = Index(lam, name ='eigenvalues')
    if flux_modes:
        columns = kmodel.reactions.keys()
    else:
        columns = [kmodel.reactants.iloc(i)[0] for i in kmodel.independent_variables_ix]

    return DataFrame(W, index=index, columns=columns)
