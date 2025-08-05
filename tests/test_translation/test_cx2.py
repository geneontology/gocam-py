import pytest
from ndex2.cx2 import CX2Network, RawCX2NetworkFactory

from gocam.translation.cx2 import model_to_cx2


def test_model_to_cx2(get_model):
    """Test the model_to_cx2 function."""
    model = get_model("input/Model-663d668500002178")
    cx2 = model_to_cx2(model)

    assert isinstance(cx2, list)

    node_aspect = next((aspect for aspect in cx2 if "nodes" in aspect), None)
    assert node_aspect is not None
    assert len(node_aspect["nodes"]) == 13, "Incorrect number of nodes in CX2"

    edge_aspect = next((aspect for aspect in cx2 if "edges" in aspect), None)
    assert edge_aspect is not None
    assert len(edge_aspect["edges"]) == 21, "Incorrect number of edges in CX2"


def test_load_cx2_to_ndex(get_model):
    """Test loading generated CX2 file by NDEx library."""
    model = get_model("input/Model-663d668500002178")
    cx2 = model_to_cx2(model)

    factory = RawCX2NetworkFactory()
    cx2_network = factory.get_cx2network(cx2)

    assert isinstance(cx2_network, CX2Network)
    assert len(cx2_network.get_nodes()) == 13, "Incorrect number of nodes in CX2"
    assert len(cx2_network.get_edges()) == 21, "Incorrect number of edges in CX2"


def test_node_type_attribute(get_model):
    """Test that the `type` attribute is correctly set for nodes."""
    model = get_model("input/Model-6606056e00002011")
    cx2 = model_to_cx2(model)

    node_aspect = next((aspect for aspect in cx2 if "nodes" in aspect), None)
    assert node_aspect is not None
    for node in node_aspect["nodes"]:
        node_attrs = node["v"]
        # If this is the expected complex node, check that the type is set to "complex"
        if node_attrs["name"] == "B cell receptor complex":
            assert node_attrs["type"] == "complex"
        # Otherwise, check that it is set to "gene" or "molecule"
        elif "type" in node_attrs:
            assert node_attrs["type"] == "gene" or node_attrs["type"] == "molecule"


def test_node_name_and_member_attributes(get_model):
    model = get_model("input/Model-6606056e00002011")
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


def test_activity_input_output_notes(get_model):
    model = get_model("input/Model-63f809ec00000701")
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
        if edge["t"] == io_node["id"] and edge["v"]["name"] == "has input"
    )
    assert input_edge is not None

    # Find the edge that has the expected source node and edge named "has output"
    output_edge = next(
        edge
        for edge in edge_aspect["edges"]
        if edge["t"] == io_node["id"] and edge["v"]["name"] == "has output"
    )
    assert output_edge is not None


def test_issue_65_protein_inputs_filtered(get_model):
    """Test that protein inputs are filtered at the CX level (issue #65)"""
    from gocam.datamodel import MoleculeAssociation

    model = get_model("input/Model-63f809ec00000701")

    # Find activities with protein inputs
    activities_with_protein_inputs = []

    # Add a protein input to an activity to test filtering
    for activity in model.activities:
        if activity.enabled_by is None:
            continue

        # Find another activity to use as input
        for other_activity in model.activities:
            if other_activity.enabled_by is None or other_activity == activity:
                continue

            # Add the other activity's enabled_by as an input to this activity
            if activity.has_input is None:
                activity.has_input = []

            activity.has_input.append(
                MoleculeAssociation(term=other_activity.enabled_by.term, evidence=[])
            )
            activities_with_protein_inputs.append(activity)
            break

        if activities_with_protein_inputs:
            break

    # Ensure we have at least one activity with a protein input for the test
    assert len(activities_with_protein_inputs) > 0

    # Convert to CX2
    cx2 = model_to_cx2(model)

    # Get edges
    edge_aspect = next((aspect for aspect in cx2 if "edges" in aspect), None)
    assert edge_aspect is not None

    # Create a list of protein terms
    protein_terms = [
        activity.enabled_by.term
        for activity in model.activities
        if activity.enabled_by is not None
    ]

    # Find all "has input" edges
    has_input_edges = [
        edge for edge in edge_aspect["edges"] if edge["v"].get("name") == "has input"
    ]

    # Get nodes
    node_aspect = next((aspect for aspect in cx2 if "nodes" in aspect), None)
    assert node_aspect is not None

    # Check if any has_input edges point to protein nodes
    for edge in has_input_edges:
        target_node = next(
            node for node in node_aspect["nodes"] if node["id"] == edge["t"]
        )

        # Verify that the target node doesn't represent a protein term
        assert target_node["v"]["represents"] not in protein_terms
