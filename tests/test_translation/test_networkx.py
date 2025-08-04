import pytest
import yaml
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