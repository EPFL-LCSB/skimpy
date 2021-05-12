Quick start
===========

In this quick start tutorial we introduce the main functions to build and manipulate SKiMpy models. 

.. code-block:: python

    import numpy as np
    from skimpy.core import *
    from skimpy.mechanisms import *

    name = 'pfk'
    metabolites = ReversibleMichaelisMenten.Reactants(substrate = 'A',
                                                       product = 'B')

    # Define reaction parameters
    parameters = ReversibleMichaelisMenten.Parameters(
        vmax_forward = 1.0,
        k_equilibrium=2.0,
        km_substrate = 10.0,
        km_product = 10.0,
        total_enzyme_concentration = 1.0,
    )

    # Create the reaction object
    pfk = Reaction(name=name,
                   mechanism = ReversibleMichaelisMenten,
                   reactants=metabolites,
                   )

    # Create a model
    this_model = KineticModel()
    
    # Add the reaction
    this_model.add_reaction(pfk)
    
    # Assing the parameters to the model 
    this_model.parametrize_by_reaction({pfk.name:parameters})
    
    # Compile the symbolic expressions of the ODE system
    this_model.compile_ode()

    # Assing initial conditions
    this_model.initial_conditions['A'] = 1.0
    this_model.initial_conditions['B'] = 1.0

    # solve the ODEs
    this_sol_qssa = this_model.solve_ode(np.linspace(0.0, 100.0, 1000), solver_type='cvode')

    this_sol_qssa.plot('output/uni_uni_base_out_qssa.html')

