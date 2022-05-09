Modal analysis
===========================

Modal analysis allows to get insight into the dynamics close the the steady state.
The modal matrix provides information on which eigenvalues contribute to the
dynamics of each concentration. Similar to MCA, Modal analysis is conducted around a steady-state
therefore modal analysis requieres to compute or import a valid set of steady-state concentrations:

.. code-block:: python

    from pytfa.io.json import load_json_model
    from skimpy.io.yaml import  load_yaml_model
    from skimpy.analysis.oracle.load_pytfa_solution import load_concentrations, load_fluxes
    from skimpy.analysis.modal import modal_matrix
    from skimpy.core.parameters import ParameterValues
    from skimpy.utils.namespace import *
    from skimpy.utils.tabdict import TabDict

    from skimpy.viz.modal import plot_modal_matrix

    # Units of the parameters are muM and hr
    CONCENTRATION_SCALING = 1e6
    TIME_SCALING = 1 # 1hr to 1min
    DENSITY = 1200 # g/L
    GDW_GWW_RATIO = 0.3 # Assumes 70% Water

    kmodel =  load_yaml_model('./../../models/varma_strain_1.yml')
    tmodel = load_json_model('./../../models/tfa_varma.json')

    # Reference steady-state data
    ref_solution = pd.read_csv('./../../data/tfa_reference_strains.csv',
                               index_col=0).loc['strain_1',:]

    ref_concentrations = load_concentrations(ref_solution, tmodel, kmodel,
                                             concentration_scaling=CONCENTRATION_SCALING)

    parameter_values = {p.symbol:p.value for p in kmodel.parameters.values()}
    parameter_values = ParameterValues(parameter_values, kmodel)



To compute the modal-matrix the model needs to have compiled Jacobian expressions, they are build by calling ``kmodel.prepare()`` and ``kmodel.compile_jacobian()``.

.. code-block:: python

    kmodel.prepare()
    kmodel.compile_jacobian(sim_type=QSSA,ncpu=8)

    M = modal_matrix(kmodel,ref_concentrations,parameter_values)

    plot_modal_matrix(M,filename='modal_matrix.html',
                      plot_width=800, plot_height=600,
                      clustered=True,
                      backend='svg',
                      )
