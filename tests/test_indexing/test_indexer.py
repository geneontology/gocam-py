"""Test for the indexer module."""

import pytest
from unittest.mock import MagicMock, patch
import networkx as nx
from linkml_runtime.loaders import yaml_loader

from gocam.datamodel import Model, QueryIndex
from gocam.indexing.indexer import Indexer
from tests import EXAMPLES_DIR, INPUT_DIR


@pytest.fixture
def example_model():
    """Load an example model for testing."""
    return yaml_loader.load(f"{EXAMPLES_DIR}/Model-663d668500002178.yaml", target_class=Model)


@pytest.fixture
def mock_go_adapter():
    """Create a mock GO adapter for testing."""
    mock_adapter = MagicMock()
    mock_adapter.ancestors.return_value = ["GO:0001", "GO:0002"]
    mock_adapter.label.return_value = "Test Term"
    return mock_adapter


def test_index_model(example_model):
    """Test that a model can be indexed."""
    # Create an indexer with a mock GO adapter
    with patch('gocam.indexing.indexer.get_adapter') as mock_get_adapter:
        mock_adapter = MagicMock()
        mock_adapter.ancestors.return_value = ["GO:0001", "GO:0002"]
        mock_adapter.label.return_value = "Test Term"
        mock_get_adapter.return_value = mock_adapter
        
        indexer = Indexer()
        
        # Test initial state
        assert example_model.query_index is None
        
        # Index the model
        indexer.index_model(example_model)
        
        # Check that query_index was created
        assert example_model.query_index is not None
        assert isinstance(example_model.query_index, QueryIndex)
        
        # Check that basic stats were calculated
        assert example_model.query_index.number_of_activities == len(example_model.activities)
        
        # Check that causal associations were counted
        assert example_model.query_index.number_of_causal_associations > 0
        
        # Check that references were flattened
        assert example_model.query_index.flattened_references is not None
        
        # Check that closures were generated for relevant terms
        assert example_model.query_index.model_activity_molecular_function_terms is not None
        assert example_model.query_index.model_activity_part_of_terms is not None
        assert example_model.query_index.model_activity_occurs_in_terms is not None
        
        # Verify reindex flag works
        # First call should not change anything since we already indexed
        indexer.index_model(example_model, reindex=False)
        # Second call with reindex=True should recreate the index
        indexer.index_model(example_model, reindex=True)


def test_model_to_digraph(example_model, mock_go_adapter):
    """Test converting a model to a directed graph."""
    with patch('gocam.indexing.indexer.get_adapter') as mock_get_adapter:
        mock_get_adapter.return_value = mock_go_adapter
        
        indexer = Indexer()
        graph = indexer.model_to_digraph(example_model)
        
        # Verify the graph structure
        assert isinstance(graph, nx.DiGraph)
        assert graph.number_of_nodes() > 0
        assert graph.number_of_edges() > 0
        
        # Verify edges correspond to causal associations
        for activity in example_model.activities:
            if activity.causal_associations:
                for ca in activity.causal_associations:
                    # Check that the edge exists from downstream activity to this activity id
                    if ca.downstream_activity:
                        assert graph.has_edge(ca.downstream_activity, activity.id)


def test_get_closures(mock_go_adapter):
    """Test getting term closures."""
    with patch('gocam.indexing.indexer.get_adapter') as mock_get_adapter:
        mock_get_adapter.return_value = mock_go_adapter
        
        indexer = Indexer()
        terms = ["GO:1234", "GO:5678"]
        
        direct, closure = indexer._get_closures(terms)
        
        # Check results
        assert len(direct) == len(terms)
        assert len(closure) == 2  # The mock returns 2 ancestors
        
        # Check structure of returned objects
        for obj in direct + closure:
            assert hasattr(obj, "id")
            assert hasattr(obj, "label")
            assert obj.label == "Test Term"  # From our mock


def test_indexer_with_empty_model():
    """Test indexer with an empty model."""
    # Create an empty model
    model = Model(id="test:empty")
    model.activities = []
    
    # Create an indexer with a mock GO adapter
    with patch('gocam.indexing.indexer.get_adapter') as mock_get_adapter:
        mock_adapter = MagicMock()
        mock_adapter.ancestors.return_value = []
        mock_adapter.label.return_value = "Test Term"
        mock_get_adapter.return_value = mock_adapter
        
        indexer = Indexer()
        indexer.index_model(model)
        
        # Verify basic properties
        assert model.query_index is not None
        assert model.query_index.number_of_activities == 0
        assert model.query_index.number_of_causal_associations == 0
        assert model.query_index.flattened_references == []
        
        # Graph should be empty
        graph = indexer.model_to_digraph(model)
        assert graph.number_of_nodes() == 0
        assert graph.number_of_edges() == 0