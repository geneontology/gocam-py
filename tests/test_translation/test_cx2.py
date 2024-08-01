import pytest
import yaml
from ndex2.cx2 import CX2Network, RawCX2NetworkFactory

from gocam.datamodel import Model
from gocam.translation.cx2 import model_to_cx2
from tests import EXAMPLES_DIR


@pytest.fixture
def model():
    model_file = EXAMPLES_DIR / "Model-663d668500002178.yaml"
    with open(model_file, "r") as f:
        deserialized = yaml.safe_load(f)
    model = Model.model_validate(deserialized)
    return model


def test_model_to_cx2(model):
    """Test the model_to_cx2 function."""
    cx2 = model_to_cx2(model)

    assert isinstance(cx2, list)

    node_aspect = next((aspect for aspect in cx2 if "nodes" in aspect), None)
    assert node_aspect is not None
    assert len(node_aspect["nodes"]) == 10, "Incorrect number of nodes in CX2"

    edge_aspect = next((aspect for aspect in cx2 if "edges" in aspect), None)
    assert edge_aspect is not None
    assert len(edge_aspect["edges"]) == 14, "Incorrect number of edges in CX2"


def test_load_cx2_to_ndex(model):
    """Test loading generated CX2 file by NDEx library."""
    cx2 = model_to_cx2(model)

    factory = RawCX2NetworkFactory()
    cx2_network = factory.get_cx2network(cx2)

    assert isinstance(cx2_network, CX2Network)
    assert len(cx2_network.get_nodes()) == 10, "Incorrect number of nodes in CX2"
    assert len(cx2_network.get_edges()) == 14, "Incorrect number of edges in CX2"
