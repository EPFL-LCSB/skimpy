from skimpy.io.yaml import export_to_yaml, load_yaml_model
from os.path import join
import pytest

#model = None#load_yaml_model()
dummy_model_path = 'test.yaml'

@pytest.mark.dependency(name='cobra_import')
def test_cobra_import():
    from skimpy.io.generate_from_cobra import FromCobra
    from cobra.io.mat import load_matlab_model

    this_cobra_model = load_matlab_model(join('..','models','toy_model.mat'), 'ToyModel_DP')

    model_gen = FromCobra()
    this_skimpy_model = model_gen.import_model(this_cobra_model)

    # Compile equations
    this_skimpy_model.compile_ode()

    # Serialize
    export_to_yaml(this_skimpy_model, dummy_model_path)


@pytest.mark.dependency(depends=['cobra_import'])
def test_import():
    model = load_yaml_model(dummy_model_path)
    model.compile_ode()
