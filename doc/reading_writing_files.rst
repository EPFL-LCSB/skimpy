Reading and writing model files
===============================
A SKiMpy model can be saved and edited using YAML, allowing users to curate large model files in a human-readable fashion.
Further, SKiMpy can generate kinetic draft models from pyTFA or cobrapy models, assigning best guess reaction mechanisms.

Function to save and load YAML model are encoded in the skimpy.io.yaml module
.. code-block:: python

    # Loading a yml model
    from skimpy.io import load_yaml_model
    kmodel = load_yaml_model('/skimpy/models/kin_varma.yml')
    
    # Exporting a SKiMpy model to yml
    from skimpy.io import export_to_yaml
    export_to_yaml(kmodel, 'test.yml')
    
    
    
Cobrapy and pyTFA models can be used in generate unparametrized draft models:
.. code-block:: python
  
  from skimpy.io.generate_from_cobra import FromCobra
  from cobra.io.mat import load_matlab_model

  # Load a cobra model
  this_cobra_model = load_matlab_model('/skimpy/models/toy_model.mat','model')

  # Generate a draft kinetic model 
  model_gen = FromCobra()
  kmodel = model_gen.import_model(this_cobra_model)

When building draft models based on cobrapy we assume full reversibility of all reactions using pyTFA models this can be specified 
by accesing the results of the Gibbs free energy vairables.

.. code-block:: python

  import pytfa
  from pytfa.io import import_matlab_model, load_thermoDB
  from pytfa.io.viz import get_reaction_data

  # Convert the cobra model to a thermodynamics model
  thermo_data = load_thermoDB('/skimpy/data/thermo_data.thermodb')
  this_pytfa_model = pytfa.ThermoModel(thermo_data, this_cobra_model)

  GLPK = 'optlang-glpk'
  this_pytfa_model.solver = GLPK

  ## TFA conversion
  this_pytfa_model.prepare()
  this_pytfa_model.convert(add_displacement=True)

  # Solve the pytfa model to find deltaGs 
  solution = this_pytfa_model.optimize()
  
  # Generate the KineticModel

  # Define the molecules that should be considered small-molecules
  # These molecules will not be accounted explicitly in the kinetic mechanism as
  # substrates and products
  small_molecules = ['h_c', 'h_e']

  model_gen = FromPyTFA(small_molecules=small_molecules)
  kmodel = model_gen.import_model(this_pytfa_model,solution.raw)
