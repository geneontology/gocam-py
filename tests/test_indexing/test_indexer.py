"""Test for the indexer module."""
import pytest
import networkx as nx
from linkml_runtime.loaders import yaml_loader

from gocam.datamodel import Model, QueryIndex
from gocam.indexing.Indexer import Indexer

from tests import EXAMPLES_DIR


@pytest.fixture
def example_model():
    """Load an example model for testing."""
    return yaml_loader.load(f"{EXAMPLES_DIR}/Model-663d668500002178.yaml", target_class=Model)


@pytest.fixture(autouse=True)
def mock_oaklib_adapters(monkeypatch):
    class MockGoAdapter:
        def ancestors(self, *args, **kwargs):
            return ["GO:0001", "GO:0002"]

        def label(self, *args, **kwargs):
            return "Test Term"

        def subset_members(self, *args, **kwargs):
            return ["GO:0001", "GO:0002", "GO:0003"]

    class MockNcbiTaxonAdapter:
        def ancestors(self, *args, **kwargs):
            return ["NCBITaxon:1", "NCBITaxon:2"]

        def label(self, *args, **kwargs):
            if args and args[0] == "NCBITaxon:9606":
                return "Human"
            return "Test Taxon"

    monkeypatch.setattr(Indexer, "go_adapter", MockGoAdapter())
    monkeypatch.setattr(Indexer, "ncbi_taxon_adapter", MockNcbiTaxonAdapter())


def test_index_model(example_model):
    """Test that a model can be indexed."""
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


def test_model_to_digraph(example_model):
    """Test converting a model to a directed graph."""
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


def test_get_closures():
    """Test getting term closures."""
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
    model = Model(id="test:empty", title="Test Empty Model")
    model.activities = []

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


def test_indexer_populates_taxon_label():
    """Test that the indexer populates taxon labels."""
    model = Model(id="test", title="Test Model", taxon="NCBITaxon:9606")
    indexer = Indexer()
    indexer.index_model(model)
    assert model.query_index.taxon_label == "Human"  # From our mock adapter


def test_indexer_gets_labels_from_model_objects():
    """Test that the indexer uses labels from the model's `objects` field when available."""
    model = Model.model_validate({
        "id": "test",
        "title": "Test Model",
        "taxon": "NCBITaxon:9606",
        "activities": [
            {
                "id": "act1",
                "enabled_by": {
                    "type": "EnabledByGeneProductAssociation",
                    "term": "UniProtKB:00001",
                }
            },
            {
                "id": "act2",
                "enabled_by": {
                    "type": "EnabledByGeneProductAssociation",
                    "term": "UniProtKB:00002",
                }
            }
        ],
        "objects": [
            {
                "id": "UniProtKB:00001",
                "label": "Test Gene"
            }
        ]
    })
    indexer = Indexer()
    indexer.index_model(model)
    assert model.query_index.model_activity_enabled_by_terms is not None
    assert len(model.query_index.model_activity_enabled_by_terms) == 2
    assert {term.label for term in model.query_index.model_activity_enabled_by_terms} == {"Test Gene", "Test Term"}


def test_indexer_adds_complex_members_to_model_activity_enabled_by_genes():
    """Test that the indexer adds complex members to model_activity_enabled_by_genes."""
    model = Model.model_validate({
        "id": "test",
        "title": "Test Model",
        "taxon": "NCBITaxon:9606",
        "activities": [
            {
                "id": "act1",
                "enabled_by": {
                    "type": "EnabledByProteinComplexAssociation",
                    "term": "GO:00001",
                    "members": ["UniProtKB:00001", "UniProtKB:00002"]
                }
            }
        ],
        "objects": [
            {
                "id": "GO:00001",
                "label": "Test Complex"
            },
            {
                "id": "UniProtKB:00001",
                "label": "Test Gene"
            },
            {
                "id": "UniProtKB:00002",
                "label": "Test Gene 2"
            }
        ]
    })
    indexer = Indexer()
    indexer.index_model(model)

    # The indexer should have unpacked the two complex members into model_activity_enabled_by_genes
    assert model.query_index.model_activity_enabled_by_genes is not None
    assert len(model.query_index.model_activity_enabled_by_genes) == 2
    assert {gene.id for gene in model.query_index.model_activity_enabled_by_genes} == {"UniProtKB:00001", "UniProtKB:00002"}
    assert {gene.label for gene in model.query_index.model_activity_enabled_by_genes} == {"Test Gene", "Test Gene 2"}

    # The original complex term should still be in model_activity_enabled_by_terms
    assert model.query_index.model_activity_enabled_by_terms is not None
    assert len(model.query_index.model_activity_enabled_by_terms) == 1
    assert model.query_index.model_activity_enabled_by_terms[0].id == "GO:00001"
    assert model.query_index.model_activity_enabled_by_terms[0].label == "Test Complex"

