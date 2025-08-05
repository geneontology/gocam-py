"""
Shared pytest fixtures for all tests.
"""
import pytest
import yaml
from pathlib import Path

from gocam.datamodel import Model
from gocam.translation.networkx.model_network_translator import ModelNetworkTranslator

# Get the tests directory path
TESTS_DIR = Path(__file__).parent
INPUT_DIR = TESTS_DIR / "input"
EXAMPLES_DIR = TESTS_DIR / "../src/data/examples"


@pytest.fixture
def get_model():
    """
    Factory fixture for loading models from YAML files.
    
    Usage:
        def test_something(get_model):
            model = get_model('/path/to/model.yaml')
    """
    def _get_model(model_path):
        with open(model_path, "r") as f:
            deserialized = yaml.safe_load(f)
        model = Model.model_validate(deserialized)
        return model
    return _get_model


@pytest.fixture
def input_model(get_model):
    """
    Factory fixture for loading models from the tests/input directory.
    
    Usage:
        def test_something(input_model):
            model = input_model('Model-63f809ec00000701')  # loads Model-63f809ec00000701.yaml
    """
    def _get_input_model(model_name):
        return get_model(INPUT_DIR / f"{model_name}.yaml")
    return _get_input_model


@pytest.fixture
def example_model(get_model):
    """
    Factory fixture for loading models from the examples directory.
    
    Usage:
        def test_something(example_model):
            model = example_model('example_model')  # loads example_model.yaml from examples/
    """
    def _get_example_model(example_name):
        return get_model(EXAMPLES_DIR / f"{example_name}.yaml")
    return _get_example_model


@pytest.fixture
def translator():
    """Create a ModelNetworkTranslator instance."""
    return ModelNetworkTranslator()