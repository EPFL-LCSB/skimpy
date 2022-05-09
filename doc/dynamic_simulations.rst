Non-linear dynamic simulations
==============================

Dynamic simulations can directly be perform by importing a model that contains
parameter values. For large-scale kinetic models though it is recommended to
parametrize the initial-conditions to a known / assumed reference steady-state
as high dimensional non-linear often exhibit multiple steady-states.


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


To run dynamic simulations the model need to contain compiled ODE expressions, by calling ``kmodel.prepare()``
and ``kmodel.compile_ode()``. To compute fluxes at specific time points it can be useful
to build a ``FluxFunction`` using the ``make_flux_fun(kmodel, QSSA)`` function.

.. code-block:: python

    kmodel.compile_ode(sim_type=QSSA,ncpu=8)
    # make function to calculate the fluxes
    flux_fun = make_flux_fun(kmodel, QSSA)

    for k in kmodel.initial_conditions:
        kmodel.initial_conditions[k] = ref_concentrations[k]

    desings = {'vmax_forward_LDH_D': 2.0,
               'vmax_forward_GAPD': 2.0}

    fluxes = []
    for p,v in desings.items():
        kmodel.parameters = parameter_values
        # Implement parameter desing
        kmodel.parameters[p].value = kmodel.parameters[p].value*v

        sol = kmodel.solve_ode(np.logspace(-9,0, 1000),
                                 solver_type='cvode')

        # Calculate fluxes:
        this_fluxes = []
        for i, concentrations in sol.concentrations.iterrows():
            t_fluxes = flux_fun(concentrations, parameters=parameter_values)
            this_fluxes.append(t_fluxes)

        # Make it a DataFrame
        this_fluxes = pd.DataFrame(this_fluxes)/ref_fluxes
        this_fluxes.index = sol.time
        fluxes.append(this_fluxes)

    ldh_fluxes = np.array([fluxes[0]['LDH_D'].values, fluxes[1]['LDH_D'].values ]).T

    timetrace_plot(sol.time,ldh_fluxes ,
                   filename='ldh_flux.html',
                   legend=['2 x [LDH]','2 x [GAPD]'],
                   backend='svg'
                   )

