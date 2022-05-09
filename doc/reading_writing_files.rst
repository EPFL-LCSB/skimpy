Reading and writing model files
===============================
A SKiMpy model can be saved and edited using YAML, allowing users to curate large model files in a human-readable fashion.
Additionally we support importing and exporting models from and to SBML format. 

Function to save and load YAML model are encoded in the skimpy.io.yaml module

.. code-block:: python

    # Loading a yml model
    from skimpy.io.yaml import load_yaml_model
    kmodel = load_yaml_model('/skimpy/models/kin_varma.yml')
    
    # Exporting a SKiMpy model to yml
    from skimpy.io.yaml import export_to_yaml
    export_to_yaml(kmodel, 'test.yml')
    
    
Function to save and load SMBL model are encoded in the skimpy.io.yaml module

.. code-block:: python

    # Exporting a SKiMpy model to SBML
    from skimpy.io.sbml import export_sbml
    export_to_yaml(kmodel, 'test.sbml')

    # Loading a SBML model
    from skimpy.io.sbml import import_sbml
    kmodel = import_sbml(kin_varma.sbml')


It should be noted here that the SMBL format does not permit to preserve the object oriented structure of the model. Thus a generic mechanism is used which cannot be exported into the prefered (Human readable) yaml format. 

