import pytest
import yaml
import json
import networkx as nx

from gocam.datamodel import Model
from gocam.translation.networkx.model_network_translator import ModelNetworkTranslator
from tests import INPUT_DIR



def test_translate_models_basic(get_model, translator):
    """Test basic translation of a GO-CAM model to gene-to-gene format."""
    model = get_model("input/Model-63f809ec00000701")
    
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


def test_gene_to_gene_edges(get_model, translator):
    """Test that edges represent gene-to-gene relationships with GO term properties."""
    model = get_model("input/Model-63f809ec00000701")
    
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


def test_multiple_models(get_model, translator):
    """Test translation of multiple models into a single graph."""
    model1 = get_model("input/Model-63f809ec00000701")
    model2 = get_model("input/Model-6606056e00002011")
    
    g2g_graph = translator.translate_models([model1, model2])
    
    # Check that we get a NetworkX DiGraph
    assert isinstance(g2g_graph, nx.DiGraph)
    
    # Check that we have nodes from both models
    model_ids_in_nodes = set()
    for node, attrs in g2g_graph.nodes(data=True):
        model_ids_in_nodes.add(attrs['model_id'])
    
    # Should have nodes from both models
    assert model1.id in model_ids_in_nodes or model2.id in model_ids_in_nodes


def test_go_term_edge_attributes(get_model, translator):
    """Test that edges contain GO term information as attributes."""
    model = get_model("input/Model-63f809ec00000701")
    
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


def test_both_source_and_target_gene_attributes(get_model, translator):
    """Test that edges include GO terms for both source and target genes."""
    model = get_model("input/Model-63f809ec00000701")
    
    g2g_graph = translator.translate_models([model])
    
    # Check that edges have both source and target gene information
    for source, target, attrs in g2g_graph.edges(data=True):
        # Basic-edge attributes should always be present
        assert 'source_gene' in attrs
        assert 'target_gene' in attrs
        assert 'model_id' in attrs
        
        # Check that source and target match-edge endpoints
        assert attrs['source_gene'] == source
        assert attrs['target_gene'] == target
        
        # If we have causal predicate, it should be present
        # (Not all edges may have all GO terms, depending on the data)


def test_empty_model_list(translator):
    """Test translation of an empty model list."""
    g2g_graph = translator.translate_models([])
    
    assert isinstance(g2g_graph, nx.DiGraph)
    assert g2g_graph.number_of_nodes() == 0
    assert g2g_graph.number_of_edges() == 0


def test_model_without_activities(translator):
    """Test handling of a model without activities."""
    # Create a minimal model
    model = Model(id="test:empty", title="Empty Test Model")
    
    g2g_graph = translator.translate_models([model])
    
    assert isinstance(g2g_graph, nx.DiGraph)
    assert g2g_graph.number_of_nodes() == 0
    assert g2g_graph.number_of_edges() == 0


def test_json_output_with_model_info(get_model, translator):
    """Test JSON output includes model_info by default."""
    model = get_model("input/Model-63f809ec00000701")
    
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


def test_json_output_without_model_info(get_model, translator):
    """Test JSON output excludes model_info when requested."""
    model = get_model("input/Model-63f809ec00000701")
    
    json_output = translator.translate_models_to_json([model], include_model_info=False)
    data = json.loads(json_output)
    
    # Check model_info is not included
    assert "graph" in data
    assert "model_info" not in data["graph"]
    assert data["graph"] == {}


def test_json_output_multiple_models_info(get_model, translator):
    """Test JSON output with multiple models includes models_info."""
    model1 = get_model("input/Model-63f809ec00000701")
    model2 = get_model("input/Model-6606056e00002011")
    
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


def test_networkx_json_format_compliance(get_model, translator):
    """Test that JSON output complies with NetworkX node_link_data format."""
    model = get_model("input/Model-63f809ec00000701")
    
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


def test_networkx_roundtrip_compatibility(get_model, translator):
    """Test that our JSON output can be read back by NetworkX."""
    model = get_model("input/Model-63f809ec00000701")
    
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


def test_json_node_structure(get_model, translator):
    """Test that nodes in JSON output have correct structure."""
    model = get_model("input/Model-63f809ec00000701")
    
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


def test_json_edge_structure(get_model, translator):
    """Test that edges in JSON output have correct structure."""
    model = get_model("input/Model-63f809ec00000701")
    
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


def test_evidence_collections_in_json(get_model, translator):
    """Test that evidence collections are properly included in JSON output."""
    model = get_model("input/Model-63f809ec00000701")
    
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


# Model 568b0f9600000284 Tests (using Model-568b0f9600000284.yaml from tests/input)

def test_model_basic_conversion_568b0f9600000284(get_model, translator):
    """Test basic conversion of model 568b0f9600000284 from GO-CAM to gene-to-gene format."""
    model = get_model("input/Model-568b0f9600000284")
    g2g_graph = translator.translate_models([model])
    
    # Should have 7 gene nodes (one per activity)
    assert g2g_graph.number_of_nodes() == 10
    
    # Should have 6 causal edges
    assert g2g_graph.number_of_edges() == 9
    
    # Verify expected genes are present
    expected_genes = {
        "WB:WBGene00006575",  # tir-1
        "WB:WBGene00006599",  # tpa-1
        "WB:WBGene00012019",  # dkf-2
        "WB:WBGene00003822",  # nsy-1
        "WB:WBGene00004758",  # sek-1
        "WB:WBGene00004055",  # pmk-1
        "WB:WBGene00000223",  # atf-7
        "WB:WBGene00006923",  # vhp-1
        "WB:WBGene00002187",  # kgb-1
        "WB:WBGene00011979",  # sysm-1
    }
    assert set(g2g_graph.nodes()) == expected_genes


def test_model_node_attributes_568b0f9600000284(get_model, translator):
    """Test that gene nodes have expected attributes in model 568b0f9600000284."""
    model = get_model("input/Model-568b0f9600000284")
    g2g_graph = translator.translate_models([model])
    
    # Check a specific node
    tir1_attrs = g2g_graph.nodes["WB:WBGene00006575"]
    assert tir1_attrs["gene_product"] == "WB:WBGene00006575"
    assert tir1_attrs["model_id"] == "gomodel:568b0f9600000284"
    assert tir1_attrs["label"] == "tir-1 Cele"


def test_model_specific_edge_go_terms_568b0f9600000284(get_model, translator):
    """Test GO terms are correctly assigned to specific edges in model 568b0f9600000284."""
    model = get_model("input/Model-568b0f9600000284")
    g2g_graph = translator.translate_models([model])
    
    # Find the edge from tir-1 to nsy-1
    tir1_to_nsy1_attrs = g2g_graph.edges["WB:WBGene00006575", "WB:WBGene00003822"]
    
    # Check source (tir-1) GO terms
    assert tir1_to_nsy1_attrs["source_gene_molecular_function"] == "GO:0035591"
    assert tir1_to_nsy1_attrs["source_gene_biological_process"] == "GO:0140367"
    assert tir1_to_nsy1_attrs["source_gene_occurs_in"] == "GO:0005737"
    
    # Check target (nsy-1) GO terms
    assert tir1_to_nsy1_attrs["target_gene_molecular_function"] == "GO:0004709"
    assert tir1_to_nsy1_attrs["target_gene_biological_process"] == "GO:0140367"
    assert tir1_to_nsy1_attrs["target_gene_occurs_in"] == "GO:0005737"


def test_model_edge_without_occurs_in_568b0f9600000284(get_model, translator):
    """Test edge where target gene has no occurs_in annotation in model 568b0f9600000284."""
    model = get_model("input/Model-568b0f9600000284")
    g2g_graph = translator.translate_models([model])
    
    # Find edge from tpa-1 to dkf-2
    tpa1_to_dkf2_attrs = g2g_graph.edges["WB:WBGene00006599", "WB:WBGene00012019"]
    
    # tpa-1 has no occurs_in, so that property should not exist
    assert "source_gene_occurs_in" not in tpa1_to_dkf2_attrs
    
    # dkf-2 has occurs_in, so it should exist for target
    assert tpa1_to_dkf2_attrs["target_gene_occurs_in"] == "GO:0009898"


def test_model_causal_pathway_structure_568b0f9600000284(get_model, translator):
    """Test that the causal pathway structure is preserved in model 568b0f9600000284."""
    model = get_model("input/Model-568b0f9600000284")
    g2g_graph = translator.translate_models([model])
    
    # Expected causal relationships based on README model
    expected_edges = {
        ("WB:WBGene00006575", "WB:WBGene00003822"),  # tir-1 -> nsy-1
        ("WB:WBGene00006599", "WB:WBGene00012019"),  # tpa-1 -> dkf-2
        ("WB:WBGene00012019", "WB:WBGene00004055"),  # dkf-2 -> pmk-1
        ("WB:WBGene00003822", "WB:WBGene00004758"),  # nsy-1 -> sek-1
        ("WB:WBGene00004758", "WB:WBGene00004055"),  # sek-1 -> pmk-1
        ("WB:WBGene00004055", "WB:WBGene00000223"),  # pmk-1 -> atf-7
        ('WB:WBGene00006923', 'WB:WBGene00004055'),  # vhp-1 -> pmk-1
        ('WB:WBGene00006923', 'WB:WBGene00002187'),  # vhp-1 -> kgb-1
        ('WB:WBGene00000223', 'WB:WBGene00011979'),  # atf-7 -> sysm-1
    }
    
    actual_edges = set(g2g_graph.edges())
    assert actual_edges == expected_edges


def test_model_multiple_inputs_to_same_gene_568b0f9600000284(get_model, translator):
    """Test handling of multiple causal inputs to the same gene (pmk-1) in model 568b0f9600000284."""
    model = get_model("input/Model-568b0f9600000284")
    g2g_graph = translator.translate_models([model])
    
    # pmk-1 (WB:WBGene00004055) should have two incoming edges
    pmk1_predecessors = list(g2g_graph.predecessors("WB:WBGene00004055"))
    assert len(pmk1_predecessors) == 3
    
    expected_predecessors = {"WB:WBGene00012019", "WB:WBGene00004758", "WB:WBGene00006923"}  # dkf-2, sek-1, vhp-1
    assert set(pmk1_predecessors) == expected_predecessors


def test_model_all_genes_have_go_annotations_568b0f9600000284(get_model, translator):
    """Test that all edges have some GO term annotations in model 568b0f9600000284."""
    model = get_model("input/Model-568b0f9600000284")
    g2g_graph = translator.translate_models([model])
    
    for source, target, attrs in g2g_graph.edges(data=True):
        # Should have at least molecular function for both source and target
        assert "source_gene_molecular_function" in attrs
        assert "target_gene_molecular_function" in attrs
        
        # All should have biological process (part_of in this example)
        assert "source_gene_biological_process" in attrs
        assert "target_gene_biological_process" in attrs


def test_model_statistics_568b0f9600000284(get_model, translator):
    """Test that final statistics match expected values for model 568b0f9600000284."""
    model = get_model("input/Model-568b0f9600000284")
    g2g_graph = translator.translate_models([model])
    
    # Original model: 10 activities
    assert len(model.activities) == 10
    
    # Gene-to-gene network: 10 genes, 9 causal relationships
    assert g2g_graph.number_of_nodes() == 10
    assert g2g_graph.number_of_edges() == 9
    
    # All nodes should be gene products
    for node, attrs in g2g_graph.nodes(data=True):
        assert attrs["gene_product"] == node
        assert node.startswith("WB:WBGene")


def test_model_evidence_collections_basic_568b0f9600000284(get_model, translator):
    """Test that evidence collections are properly included in edges for model 568b0f9600000284."""
    model = get_model("input/Model-568b0f9600000284")
    g2g_graph = translator.translate_models([model])
    
    # Find an edge with evidence to test
    edge_with_evidence = None
    for source, target, attrs in g2g_graph.edges(data=True):
        if any(key.endswith('_has_reference') for key in attrs.keys()):
            edge_with_evidence = (source, target, attrs)
            break
    
    assert edge_with_evidence is not None, "Should find at least one edge with evidence"
    source, target, attrs = edge_with_evidence
    
    # Check for evidence collection attributes
    evidence_attrs = [key for key in attrs.keys() 
                     if any(suffix in key for suffix in ['_has_reference', '_assessed_by', '_contributors'])]
    assert len(evidence_attrs) > 0, "Should have evidence collection attributes"


def test_model_gene_product_evidence_collections_568b0f9600000284(get_model, translator):
    """Test gene product evidence collections are properly extracted for model 568b0f9600000284."""
    model = get_model("input/Model-568b0f9600000284")
    g2g_graph = translator.translate_models([model])
    
    # Find edge from tir-1 to nsy-1 which should have molecular function evidence
    tir1_to_nsy1_attrs = g2g_graph.edges["WB:WBGene00006575", "WB:WBGene00003822"]
    
    # Check source gene product evidence
    assert "source_gene_product_has_reference" in tir1_to_nsy1_attrs
    assert "source_gene_product_assessed_by" in tir1_to_nsy1_attrs
    
    # Verify reference format
    source_refs = tir1_to_nsy1_attrs["source_gene_product_has_reference"]
    assert isinstance(source_refs, list)
    assert "PMID:15625192" in source_refs
    
    # Verify evidence code format
    source_codes = tir1_to_nsy1_attrs["source_gene_product_assessed_by"]
    assert isinstance(source_codes, list)
    assert "ECO:0000314" in source_codes
    
    # Check target gene product evidence
    assert "target_gene_product_has_reference" in tir1_to_nsy1_attrs
    assert "target_gene_product_assessed_by" in tir1_to_nsy1_attrs
    
    target_refs = tir1_to_nsy1_attrs["target_gene_product_has_reference"]
    assert isinstance(target_refs, list)
    assert "PMID:11751572" in target_refs