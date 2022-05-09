Metabolic control analysis
===========================


Further we need to compute the steady state concentrations and fluxes for the steady-state
we aim to perform the MCA for. Here we use the concentrations and fluxes
from we used for the respective ORACLE parametrization.

.. code-block:: python

    from pytfa.io.json import load_json_model
    from skimpy.io.yaml import  load_yaml_model
    from skimpy.analysis.oracle.load_pytfa_solution import load_concentrations, load_fluxes
    from skimpy.core.parameters import ParameterValues
    from skimpy.utils.namespace import *

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
    ref_fluxes = load_fluxes(ref_solution, tmodel, kmodel,
                                   density=DENSITY,
                                   ratio_gdw_gww=GDW_GWW_RATIO,
                                   concentration_scaling=CONCENTRATION_SCALING,
                                   time_scaling=TIME_SCALING)

    # Extract the current parameter set from the model
    parameter_values = {p.symbol:p.value for p in kmodel.parameters.values()}
    parameter_values = ParameterValues(parameter_values, kmodel)


With a set of parameters, fluxes and concentration we can then calculate the control coefficients.
To perform metabolic control analysis the control-coeffcient expressions
need to be compiled by calling ``kmodel.prepare()`` and ``kmodel.compile_mca(parameter_list=parameter_list)``,
where ``parameter_list`` is a TabDict contaning the symbols of the parameters we
want to the control-coeffcients for. Here we calculate the control-coeffcients with respect
to all maximal enzyme rates ``V_max``.


.. code-block:: python

    from skimpy.utils.tabdict import TabDict
    from skimpy.viz.controll_coefficients import plot_control_coefficients


    # Compile mca with parameter elasticities with respect to Vmaxes
    parameter_list = TabDict([(k, p.symbol) for k, p in
                              kmodel.parameters.items() if
                              p.name.startswith('vmax_forward')])

    kmodel.compile_mca(sim_type=QSSA,ncpu=8, parameter_list=parameter_list)

    flux_control_coeff = kmodel.flux_control_fun(ref_fluxes,
                                                 ref_concentrations,
                                                 [parameter_values, ])

    lac_control_coeff = flux_control_coeff.slice_by('sample',0).loc['LDH_D', :]

    lac_control_coeff.index = [v.replace('vmax_forward_','')
                                for v in lac_control_coeff.index ]
    plot_control_coefficients(lac_control_coeff,
                              filename='lac_control_coeff.html',
                              backend='svg',
                              )


