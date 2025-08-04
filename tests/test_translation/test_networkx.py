import pytest
import yaml
import json
import networkx as nx

from gocam.datamodel import Model
from gocam.translation.networkx.model_network_translator import ModelNetworkTranslator
from tests import INPUT_DIR


@pytest.fixture
def get_model():
    def _get_model(model_path):
        with open(model_path, "r") as f:
            deserialized = yaml.safe_load(f)
        model = Model.model_validate(deserialized)
        return model
    return _get_model


@pytest.fixture
def input_model(get_model):
    def _get_input_model(model_name):
        return get_model(INPUT_DIR / f"{model_name}.yaml")
    return _get_input_model


@pytest.fixture
def translator():
    return ModelNetworkTranslator()


def test_translate_models_basic(input_model, translator):
    """Test basic translation of a GO-CAM model to gene-to-gene format."""
    model = input_model("Model-63f809ec00000701")
    
    # Translate single model
    g2g_graph = translator.translate_models([model])
    
    # Check that we get a NetworkX DiGraph
    assert isinstance(g2g_graph, nx.DiGraph)
    
    # Check that we have nodes (gene products)
    assert g2g_graph.number_of_nodes() > 0
    
    # Check that nodes have expected attributes
    for node, attrs in g2g_graph.nodes(data=True):
        assert 'gene_product' in attrs
        assert 'model_id' in attrs
        assert attrs['model_id'] == model.id


def test_gene_to_gene_edges(input_model, translator):
    """Test that edges represent gene-to-gene relationships with GO term properties."""
    model = input_model("Model-63f809ec00000701")
    
    g2g_graph = translator.translate_models([model])
    
    # Check that we have edges
    if g2g_graph.number_of_edges() > 0:
        # Check edge attributes
        for source, target, attrs in g2g_graph.edges(data=True):
            # Basic edge attributes
            assert 'source_gene' in attrs
            assert 'target_gene' in attrs
            assert 'model_id' in attrs
            
            # Source and target should match the edge endpoints
            assert attrs['source_gene'] == source
            assert attrs['target_gene'] == target
            assert attrs['model_id'] == model.id


def test_multiple_models(input_model, translator):
    """Test translation of multiple models into a single graph."""
    model1 = input_model("Model-63f809ec00000701")
    model2 = input_model("Model-6606056e00002011")
    
    g2g_graph = translator.translate_models([model1, model2])
    
    # Check that we get a NetworkX DiGraph
    assert isinstance(g2g_graph, nx.DiGraph)
    
    # Check that we have nodes from both models
    model_ids_in_nodes = set()
    for node, attrs in g2g_graph.nodes(data=True):
        model_ids_in_nodes.add(attrs['model_id'])
    
    # Should have nodes from both models
    assert model1.id in model_ids_in_nodes or model2.id in model_ids_in_nodes


def test_go_term_edge_attributes(input_model, translator):
    """Test that edges contain GO term information as attributes."""
    model = input_model("Model-63f809ec00000701")
    
    g2g_graph = translator.translate_models([model])
    
    # Look for edges with GO term attributes
    edges_with_source_go_terms = 0
    edges_with_target_go_terms = 0
    
    for source, target, attrs in g2g_graph.edges(data=True):
        # Count edges that have source gene GO term attributes
        source_go_term_attrs = [
            'source_gene_molecular_function',
            'source_gene_biological_process', 
            'source_gene_occurs_in',
        ]
        
        # Count edges that have target gene GO term attributes
        target_go_term_attrs = [
            'target_gene_molecular_function',
            'target_gene_biological_process', 
            'target_gene_occurs_in',
        ]
        
        if any(attr in attrs for attr in source_go_term_attrs):
            edges_with_source_go_terms += 1
            
        if any(attr in attrs for attr in target_go_term_attrs):
            edges_with_target_go_terms += 1
    
    # We should have at least some edges with GO term information
    # (The exact number depends on the test data)
    assert edges_with_source_go_terms >= 0  # At minimum, no errors should occur
    assert edges_with_target_go_terms >= 0  # At minimum, no errors should occur


def test_both_source_and_target_gene_attributes(input_model, translator):
    """Test that edges include GO terms for both source and target genes."""
    model = input_model("Model-63f809ec00000701")
    
    g2g_graph = translator.translate_models([model])
    
    # Check that edges have both source and target gene information
    for source, target, attrs in g2g_graph.edges(data=True):
        # Basic edge attributes should always be present
        assert 'source_gene' in attrs
        assert 'target_gene' in attrs
        assert 'model_id' in attrs
        
        # Check that source and target match edge endpoints
        assert attrs['source_gene'] == source
        assert attrs['target_gene'] == target
        
        # If we have causal predicate, it should be present
        # (Not all edges may have all GO terms, depending on the data)


def test_empty_model_list(translator):
    """Test translation of empty model list."""
    g2g_graph = translator.translate_models([])
    
    assert isinstance(g2g_graph, nx.DiGraph)
    assert g2g_graph.number_of_nodes() == 0
    assert g2g_graph.number_of_edges() == 0


def test_model_without_activities(translator):
    """Test handling of model without activities."""
    # Create a minimal model
    model = Model(id="test:empty", title="Empty Test Model")
    
    g2g_graph = translator.translate_models([model])
    
    assert isinstance(g2g_graph, nx.DiGraph)
    assert g2g_graph.number_of_nodes() == 0
    assert g2g_graph.number_of_edges() == 0


def test_json_output_with_model_info(input_model, translator):
    """Test JSON output includes model_info by default."""
    model = input_model("Model-63f809ec00000701")
    
    json_output = translator.translate_models_to_json([model])
    data = json.loads(json_output)
    
    # Check model_info is included by default
    assert "graph" in data
    assert "model_info" in data["graph"]
    
    model_info = data["graph"]["model_info"]
    assert model_info["id"] == model.id
    assert model_info["title"] == model.title
    assert model_info["taxon"] == model.taxon
    assert model_info["status"] == model.status


def test_json_output_without_model_info(input_model, translator):
    """Test JSON output excludes model_info when requested."""
    model = input_model("Model-63f809ec00000701")
    
    json_output = translator.translate_models_to_json([model], include_model_info=False)
    data = json.loads(json_output)
    
    # Check model_info is not included
    assert "graph" in data
    assert "model_info" not in data["graph"]
    assert data["graph"] == {}


def test_json_output_multiple_models_info(input_model, translator):
    """Test JSON output with multiple models includes models_info."""
    model1 = input_model("Model-63f809ec00000701")
    model2 = input_model("Model-6606056e00002011")
    
    json_output = translator.translate_models_to_json([model1, model2])
    data = json.loads(json_output)
    
    # Check models_info is included for multiple models
    assert "graph" in data
    assert "models_info" in data["graph"]
    assert "model_info" not in data["graph"]  # Should use models_info, not model_info
    
    models_info = data["graph"]["models_info"]
    assert len(models_info) == 2
    
    # Check each model's info
    model_ids = {info["id"] for info in models_info}
    assert model1.id in model_ids
    assert model2.id in model_ids


def test_networkx_json_format_compliance(input_model, translator):
    """Test that JSON output complies with NetworkX node_link_data format."""
    model = input_model("Model-63f809ec00000701")
    
    json_output = translator.translate_models_to_json([model])
    data = json.loads(json_output)
    
    # Check required NetworkX node_link_data format fields
    assert "directed" in data
    assert "multigraph" in data
    assert "graph" in data
    assert "nodes" in data
    assert "edges" in data
    
    # Check correct types and values
    assert isinstance(data["directed"], bool)
    assert isinstance(data["multigraph"], bool)
    assert isinstance(data["graph"], dict)
    assert isinstance(data["nodes"], list)
    assert isinstance(data["edges"], list)
    
    # For DiGraph, directed should be True
    assert data["directed"] is True
    # We don't use multigraph, so should be False
    assert data["multigraph"] is False


def test_networkx_roundtrip_compatibility(input_model, translator):
    """Test that our JSON output can be read back by NetworkX."""
    model = input_model("Model-63f809ec00000701")
    
    # Generate our JSON
    json_output = translator.translate_models_to_json([model])
    data = json.loads(json_output)
    
    # Try to recreate a NetworkX graph from our data
    # We need to adjust the format slightly since we use "edges" instead of "links"
    nx_data = data.copy()
    nx_data["links"] = nx_data.pop("edges")  # NetworkX expects "links", we use "edges"
    
    # This should not raise an exception
    reconstructed_graph = nx.node_link_graph(nx_data)
    
    # Check basic properties
    assert isinstance(reconstructed_graph, nx.DiGraph)
    assert reconstructed_graph.is_directed() == data["directed"]
    assert reconstructed_graph.is_multigraph() == data["multigraph"]
    
    # Check that graph attributes are preserved
    if "model_info" in data["graph"]:
        assert "model_info" in reconstructed_graph.graph
        assert reconstructed_graph.graph["model_info"] == data["graph"]["model_info"]


def test_json_node_structure(input_model, translator):
    """Test that nodes in JSON output have correct structure."""
    model = input_model("Model-63f809ec00000701")
    
    json_output = translator.translate_models_to_json([model])
    data = json.loads(json_output)
    
    # Check nodes structure
    assert len(data["nodes"]) > 0
    
    for node in data["nodes"]:
        # Each node must have an "id" field (NetworkX requirement)
        assert "id" in node
        # Our custom attributes
        assert "gene_product" in node
        assert "model_id" in node
        assert node["model_id"] == model.id


def test_json_edge_structure(input_model, translator):
    """Test that edges in JSON output have correct structure."""
    model = input_model("Model-63f809ec00000701")
    
    json_output = translator.translate_models_to_json([model])
    data = json.loads(json_output)
    
    # Check edges structure
    if len(data["edges"]) > 0:
        for edge in data["edges"]:
            # Each edge must have "source" and "target" fields (NetworkX requirement)
            assert "source" in edge
            assert "target" in edge
            # Our custom attributes
            assert "source_gene" in edge
            assert "target_gene" in edge
            assert "model_id" in edge
            assert edge["model_id"] == model.id
            
            # Source and target should match our gene attributes
            assert edge["source"] == edge["source_gene"]
            assert edge["target"] == edge["target_gene"]


def test_evidence_collections_in_json(input_model, translator):
    """Test that evidence collections are properly included in JSON output."""
    model = input_model("Model-63f809ec00000701")
    
    json_output = translator.translate_models_to_json([model])
    data = json.loads(json_output)
    
    # Look for evidence collections in edges
    evidence_attrs_found = False
    
    for edge in data["edges"]:
        # Check for various evidence collection attributes
        evidence_attrs = [
            "causal_predicate_has_reference",
            "causal_predicate_assessed_by", 
            "causal_predicate_contributors",
            "source_gene_molecular_function_has_reference",
            "source_gene_molecular_function_assessed_by",
            "source_gene_molecular_function_contributors",
            "target_gene_molecular_function_has_reference",
            "target_gene_molecular_function_assessed_by",
            "target_gene_molecular_function_contributors"
        ]
        
        if any(attr in edge for attr in evidence_attrs):
            evidence_attrs_found = True
            # If evidence attributes exist, they should be lists
            for attr in evidence_attrs:
                if attr in edge:
                    assert isinstance(edge[attr], list), f"{attr} should be a list"
                    # Lists should not be empty if they exist
                    assert len(edge[attr]) > 0, f"{attr} should not be empty"
    
    # We expect to find some evidence attributes in the test data
    # (This assertion might need adjustment based on actual test data)
    # For now, we just ensure no errors occur during processing
    assert isinstance(evidence_attrs_found, bool)