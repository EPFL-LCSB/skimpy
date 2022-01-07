Quick start
===========

In this quick start tutorial we introduce the main functions ODE to build and manipulate SKiMpy models using oversimplfied model containing one reaction. 

All main modelling functionality is imported from submodules of 'skimpy.core' and the libary of mechanisms is imported from 'skimpy.mechanisms'. The model is initialized by the 'KineticModel' constructor. 

.. code-block:: python

    # Load core functionalities + numpy (because useful)
    import numpy as np
    from skimpy.core import *
    from skimpy.mechanisms import *
    
     # Create a model
    this_model = KineticModel()

Reaction objects are constricted from a name, a mechanism object and a reactants set inherited from the respective mechanism object to ensure the consitency of e.g. order kinetic mechanisms.

.. code-block:: python

    # Name the reaction
    name = 'pfk'
    
    # Define reaction reactants sets
    metabolites = ReversibleMichaelisMenten.Reactants(substrate = 'A',
                                                       product = 'B')
    
    # Create the reaction object
    pfk = Reaction(name=name,
                   mechanism = ReversibleMichaelisMenten,
                   reactants=metabolites,
                   )
    
    # Add the reaction
    this_model.add_reaction(pfk)
    
Each mechanism has respective Parameter "ItemSets' that defines the respective parameters for the mechanism. Mechanism specific parameters set can be created using the inherited 
Mechanism.Parameters constructor, and can then be assigned using the model.parametrize_by_reaction() method: 

.. code-block:: python

    # Define reaction parameters
    parameters = ReversibleMichaelisMenten.Parameters(
        vmax_forward = 1.0,
        k_equilibrium=2.0,
        km_substrate = 10.0,
        km_product = 10.0,
        total_enzyme_concentration = 1.0,
    )
    
    # Assing the parameters to the model 
    this_model.parametrize_by_reaction({pfk.name:parameters})
 
Alternatively parameters can be directly assinged using the model.parameter poperty:

.. code-block:: python

    # Direct assginment of parameter values
    this_model.parameters.vmax_forward_pfk.value = 2.0

In this way simple models can be constructed. To integrate the ordinary differntial equations we first compile the respective functions using CYTHON by calling the 
KineticModel.compile_ode() method. Then the initial value problem can be initialized and the model can be solved using the KineticModel.solve_ode() method.
Additionally, the solution object has standardized ploting methods implemented (Note that an output folde needs to be generated manually).

.. code-block:: python

    # Compile the symbolic expressions of the ODE system
    this_model.compile_ode()

    # Assing initial conditions
    this_model.initial_conditions['A'] = 1.0
    this_model.initial_conditions['B'] = 1.0

    # solve the ODEs
    this_sol_qssa = this_model.solve_ode(np.linspace(0.0, 100.0, 1000), solver_type='cvode')

    this_sol_qssa.plot('output/uni_uni_base_out_qssa.html')

