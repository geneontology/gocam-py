import pytest
import yaml
from ndex2.cx2 import CX2Network, RawCX2NetworkFactory

from gocam.datamodel import Model
from gocam.translation.cx2 import model_to_cx2
from tests import EXAMPLES_DIR, INPUT_DIR


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


def test_model_to_cx2(example_model):
    """Test the model_to_cx2 function."""
    model = example_model("Model-663d668500002178")
    cx2 = model_to_cx2(model)

    assert isinstance(cx2, list)

    node_aspect = next((aspect for aspect in cx2 if "nodes" in aspect), None)
    assert node_aspect is not None
    assert len(node_aspect["nodes"]) == 10, "Incorrect number of nodes in CX2"

    edge_aspect = next((aspect for aspect in cx2 if "edges" in aspect), None)
    assert edge_aspect is not None
    assert len(edge_aspect["edges"]) == 14, "Incorrect number of edges in CX2"


def test_load_cx2_to_ndex(example_model):
    """Test loading generated CX2 file by NDEx library."""
    model = example_model("Model-663d668500002178")
    cx2 = model_to_cx2(model)

    factory = RawCX2NetworkFactory()
    cx2_network = factory.get_cx2network(cx2)

    assert isinstance(cx2_network, CX2Network)
    assert len(cx2_network.get_nodes()) == 10, "Incorrect number of nodes in CX2"
    assert len(cx2_network.get_edges()) == 14, "Incorrect number of edges in CX2"


def test_node_type_attribute(input_model):
    """Test that the `type` attribute is correctly set for nodes."""
    model = input_model("Model-6606056e00002011")
    cx2 = model_to_cx2(model)

    node_aspect = next((aspect for aspect in cx2 if "nodes" in aspect), None)
    assert node_aspect is not None
    for node in node_aspect["nodes"]:
        node_attrs = node["v"]
        # If this is the expected complex node, check that the type is set to "complex"
        if node_attrs["name"] == "B cell receptor complex":
            assert node_attrs["type"] == "complex"
        # Otherwise, if the node has a type attribute (nodes created by input/output
        # associations won't have that attribute), check that it is set to "gene"
        elif "type" in node_attrs:
            assert node_attrs["type"] == "gene"


def test_node_name_and_member_attributes(input_model):
    model = input_model("Model-6606056e00002011")
    cx2 = model_to_cx2(model)

    node_aspect = next((aspect for aspect in cx2 if "nodes" in aspect), None)
    assert node_aspect is not None
    for node in node_aspect["nodes"]:
        node_attrs = node["v"]
        if node_attrs["name"] == "B cell receptor complex":
            assert "member" in node_attrs
            assert len(node_attrs["member"]) == 2
            assert all("Hsap" not in member for member in node_attrs["member"])
        else:
            assert "member" not in node_attrs
            assert "Hsap" not in node_attrs["name"]


def test_activity_input_output_notes(input_model):
    model = input_model("Model-63f809ec00000701")
    cx2 = model_to_cx2(model)

    node_aspect = next((aspect for aspect in cx2 if "nodes" in aspect), None)
    assert node_aspect is not None
    # Find the node that should be the source of both a "has input" and "has output" edge
    io_node = next(
        (
            node
            for node in node_aspect["nodes"]
            if node["v"]["name"] == "tRNA precursor"
        ),
        None,
    )
    assert io_node is not None

    edge_aspect = next((aspect for aspect in cx2 if "edges" in aspect), None)
    assert edge_aspect is not None

    # Find the edge that has the expected source node and edge named "has input"
    input_edge = next(
        edge
        for edge in edge_aspect["edges"]
        if edge["s"] == io_node["id"] and edge["v"]["name"] == "has input"
    )
    assert input_edge is not None

    # Find the edge that has the expected source node and edge named "has output"
    output_edge = next(
        edge
        for edge in edge_aspect["edges"]
        if edge["s"] == io_node["id"] and edge["v"]["name"] == "has output"
    )
    assert output_edge is not None
