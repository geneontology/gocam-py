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


@pytest.fixture
def get_model():
    """
    Factory fixture for loading models from YAML files.
    
    Usage:
        def test_something(get_model):
            model = get_model('input/Model-63f809ec00000701')  # loads from tests/input/
    """
    def _get_model(model_path):
        if model_path.startswith('input/'):
            full_path = INPUT_DIR / f"{model_path[6:]}.yaml"
        else:
            raise ValueError(f"Model path must start with 'input/', got: {model_path}")
        
        with open(full_path, "r") as f:
            deserialized = yaml.safe_load(f)
        model = Model.model_validate(deserialized)
        return model
    return _get_model


@pytest.fixture
def translator():
    """Create a ModelNetworkTranslator instance."""
    return ModelNetworkTranslator()