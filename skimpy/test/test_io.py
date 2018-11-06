from .io.yaml import export_to_yaml, load_yaml_model


model = None#load_yaml_model()

def test_export_yaml:
    export_to_yaml(model)