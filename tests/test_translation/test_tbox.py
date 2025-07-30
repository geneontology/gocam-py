import pytest
import yaml

from gocam.datamodel import Model
from gocam.translation.tbox_translator import TBoxTranslator
from tests import EXAMPLES_DIR, INPUT_DIR, OUTPUT_DIR


# TODO: DRY
@pytest.fixture
def get_model():
    def _get_model(model_path):
        with open(model_path, "r") as f:
            deserialized = yaml.safe_load(f)
        model = Model.model_validate(deserialized)
        return model

    return _get_model


@pytest.fixture
def example_model(get_model):
    def _get_example_model(example_name):
        return get_model(EXAMPLES_DIR / f"{example_name}.yaml")

    return _get_example_model


@pytest.fixture
def input_model(get_model):
    def _get_input_model(model_name):
        return get_model(INPUT_DIR / f"{model_name}.yaml")

    return _get_input_model

def test_model_to_tbox(example_model):
    """Test the model_to_cx2 function."""
    tbox_translator = TBoxTranslator()
    model_acc = "Model-663d668500002178"
    model = example_model(model_acc)
    tbox_translator.load_models([model])
    ont = tbox_translator.ontology
    axioms = ont.get_axioms()
    for a in axioms:
        print(a)
    OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
    ont.save_to_file(str(OUTPUT_DIR / f"{model_acc}.ofn"))
