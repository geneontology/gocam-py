#!/usr/bin/env python3
"""
Test gene-to-gene serialization using the README GO-CAM example.
"""
import json
import yaml
import pytest
from gocam.datamodel import Model
from gocam.translation.networkx.model_network_translator import ModelNetworkTranslator


# GO-CAM model from README.md (simplified version for testing)
README_GOCAM_YAML = """
id: gomodel:568b0f9600000284
title: Antibacterial innate immune response in the intestine via MAPK cascade (C. elegans)
taxon: NCBITaxon:6239
status: production
comments:
- 'Automated change 2023-03-16: RO:0002212 replaced by RO:0002630'
- 'Automated change 2023-03-16: RO:0002213 replaced by RO:0002629'
activities:
- id: gomodel:568b0f9600000284/57ec3a7e00000079
  enabled_by:
    term: WB:WBGene00006575
  molecular_function:
    evidence:
    - term: ECO:0000314
      reference: PMID:15625192
      provenances:
      - contributor: 
        - https://orcid.org/0000-0002-1706-4196
        date: '2019-09-23'
    term: GO:0035591
  occurs_in:
    evidence:
    - term: ECO:0000314
      reference: PMID:15625192
      provenances:
      - contributor: 
        - https://orcid.org/0000-0002-3013-9906
        date: '2019-06-28'
    term: GO:0005737
  part_of:
    evidence:
    - term: ECO:0000315
      reference: PMID:19837372
      provenances:
      - contributor: 
        - https://orcid.org/0000-0002-1706-4196
        date: '2021-07-08'
    term: GO:0140367
  causal_associations:
  - evidence:
    - term: ECO:0000315
      reference: PMID:15123841
      provenances:
      - contributor: 
        - https://orcid.org/0000-0002-3013-9906
        date: '2019-05-31'
    predicate: RO:0002629
    downstream_activity: gomodel:568b0f9600000284/57ec3a7e00000109
- id: gomodel:568b0f9600000284/5b528b1100002286
  enabled_by:
    term: WB:WBGene00006599
  molecular_function:
    evidence:
    - term: ECO:0000501
      reference: GO_REF:0000037
    term: GO:0004674
  part_of:
    evidence:
    - term: ECO:0000315
      reference: PMID:22470487
    term: GO:0002225
  causal_associations:
  - evidence:
    - term: ECO:0000316
      reference: PMID:19371715
    predicate: RO:0002629
    downstream_activity: gomodel:568b0f9600000284/5745387b00000588
- id: gomodel:568b0f9600000284/5745387b00000588
  enabled_by:
    term: WB:WBGene00012019
  molecular_function:
    evidence:
    - term: ECO:0000314
      reference: PMID:17728253
    term: GO:0004674
  occurs_in:
    evidence:
    - term: ECO:0000314
      reference: PMID:17728253
    term: GO:0009898
  part_of:
    evidence:
    - term: ECO:0000315
      reference: PMID:19371715
    term: GO:0140367
  causal_associations:
  - evidence:
    - term: ECO:0000315
      reference: PMID:19371715
    predicate: RO:0002629
    downstream_activity: gomodel:568b0f9600000284/568b0f9600000285
- id: gomodel:568b0f9600000284/57ec3a7e00000109
  enabled_by:
    term: WB:WBGene00003822
  molecular_function:
    evidence:
    - term: ECO:0000314
      reference: PMID:11751572
    term: GO:0004709
  occurs_in:
    evidence:
    - term: ECO:0000314
      reference: PMID:15625192
    term: GO:0005737
  part_of:
    evidence:
    - term: ECO:0000315
      reference: PMID:12142542
    term: GO:0140367
  causal_associations:
  - evidence:
    - term: ECO:0000315
      reference: PMID:12142542
    predicate: RO:0002629
    downstream_activity: gomodel:568b0f9600000284/57ec3a7e00000119
- id: gomodel:568b0f9600000284/57ec3a7e00000119
  enabled_by:
    term: WB:WBGene00004758
  molecular_function:
    evidence:
    - term: ECO:0000314
      reference: PMID:11751572
    term: GO:0004708
  occurs_in:
    evidence:
    - term: ECO:0000318
      reference: PMID:21873635
    term: GO:0005737
  part_of:
    evidence:
    - term: ECO:0000315
      reference: PMID:12142542
    term: GO:0140367
  causal_associations:
  - evidence:
    - term: ECO:0000315
      reference: PMID:12142542
    predicate: RO:0002629
    downstream_activity: gomodel:568b0f9600000284/568b0f9600000285
- id: gomodel:568b0f9600000284/568b0f9600000285
  enabled_by:
    term: WB:WBGene00004055
  molecular_function:
    evidence:
    - term: ECO:0000314
      reference: PMID:20369020
    term: GO:0004707
  occurs_in:
    evidence:
    - term: ECO:0000314
      reference: PMID:20133945
    term: GO:0005829
  part_of:
    evidence:
    - term: ECO:0000315
      reference: PMID:12142542
    term: GO:0140367
  causal_associations:
  - evidence:
    - term: ECO:0000314
      reference: PMID:20369020
    predicate: RO:0002629
    downstream_activity: gomodel:568b0f9600000284/568b0f9600000287
- id: gomodel:568b0f9600000284/568b0f9600000287
  enabled_by:
    term: WB:WBGene00000223
  molecular_function:
    evidence:
    - term: ECO:0000250
      reference: PMID:20369020
    term: GO:0000981
  occurs_in:
    evidence:
    - term: ECO:0000314
      reference: PMID:20369020
    term: GO:0005634
  part_of:
    evidence:
    - term: ECO:0000315
      reference: PMID:20369020
    term: GO:0140367
objects:
- id: WB:WBGene00006575
  label: tir-1 Cele
- id: WB:WBGene00006599
  label: tpa-1 Cele
- id: WB:WBGene00012019
  label: dkf-2 Cele
- id: WB:WBGene00003822
  label: nsy-1 Cele
- id: WB:WBGene00004758
  label: sek-1 Cele
- id: WB:WBGene00004055
  label: pmk-1 Cele
- id: WB:WBGene00000223
  label: atf-7 Cele
- id: GO:0035591
  label: signaling adaptor activity
- id: GO:0004674
  label: protein serine/threonine kinase activity
- id: GO:0004709
  label: MAP kinase kinase kinase activity
- id: GO:0004708
  label: MAP kinase kinase activity
- id: GO:0004707
  label: MAP kinase activity
- id: GO:0000981
  label: DNA-binding transcription factor activity, RNA polymerase II-specific
- id: GO:0005737
  label: cytoplasm
- id: GO:0009898
  label: cytoplasmic side of plasma membrane
- id: GO:0005829
  label: cytosol
- id: GO:0005634
  label: nucleus
- id: GO:0140367
  label: antibacterial innate immune response
- id: GO:0002225
  label: positive regulation of antimicrobial peptide production
- id: RO:0002629
  label: directly positively regulates
"""


@pytest.fixture
def readme_model():
    """Create a GO-CAM model from the README example."""
    gocam_data = yaml.safe_load(README_GOCAM_YAML)
    return Model.model_validate(gocam_data)


@pytest.fixture
def translator():
    """Create a ModelNetworkTranslator instance."""
    return ModelNetworkTranslator()


class TestG2GSerialization:
    """Test suite for gene-to-gene serialization using the README example."""
    
    def test_basic_g2g_conversion(self, readme_model, translator):
        """Test basic conversion from GO-CAM to gene-to-gene format."""
        g2g_graph = translator.translate_models([readme_model])
        
        # Should have 7 gene nodes (one per activity)
        assert g2g_graph.number_of_nodes() == 7
        
        # Should have 6 causal edges
        assert g2g_graph.number_of_edges() == 6
        
        # Verify expected genes are present
        expected_genes = {
            "WB:WBGene00006575",  # tir-1
            "WB:WBGene00006599",  # tpa-1
            "WB:WBGene00012019",  # dkf-2
            "WB:WBGene00003822",  # nsy-1
            "WB:WBGene00004758",  # sek-1
            "WB:WBGene00004055",  # pmk-1
            "WB:WBGene00000223",  # atf-7
        }
        assert set(g2g_graph.nodes()) == expected_genes
    
    def test_node_attributes(self, readme_model, translator):
        """Test that gene nodes have expected attributes."""
        g2g_graph = translator.translate_models([readme_model])
        
        # Check a specific node
        tir1_attrs = g2g_graph.nodes["WB:WBGene00006575"]
        assert tir1_attrs["gene_product"] == "WB:WBGene00006575"
        assert tir1_attrs["model_id"] == "gomodel:568b0f9600000284"
        assert tir1_attrs["label"] == "tir-1 Cele"
    
    def test_edge_attributes_basic(self, readme_model, translator):
        """Test that edges have basic required attributes."""
        g2g_graph = translator.translate_models([readme_model])
        
        # Get any edge to test basic structure
        edge = list(g2g_graph.edges(data=True))[0]
        source, target, attrs = edge
        
        # Basic required attributes
        assert "source_gene" in attrs
        assert "target_gene" in attrs
        assert "model_id" in attrs
        assert "causal_predicate" in attrs
        
        # Should be RO:0002629 for all edges in this example
        assert attrs["causal_predicate"] == "RO:0002629"
    
    def test_specific_edge_go_terms(self, readme_model, translator):
        """Test GO terms are correctly assigned to specific edges."""
        g2g_graph = translator.translate_models([readme_model])
        
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
    
    def test_edge_without_occurs_in(self, readme_model, translator):
        """Test edge where target gene has no occurs_in annotation."""
        g2g_graph = translator.translate_models([readme_model])
        
        # Find edge from tpa-1 to dkf-2
        tpa1_to_dkf2_attrs = g2g_graph.edges["WB:WBGene00006599", "WB:WBGene00012019"]
        
        # tpa-1 has no occurs_in, so that property should not exist
        assert "source_gene_occurs_in" not in tpa1_to_dkf2_attrs
        
        # dkf-2 has occurs_in, so it should exist for target
        assert tpa1_to_dkf2_attrs["target_gene_occurs_in"] == "GO:0009898"
    
    def test_json_serialization(self, readme_model, translator):
        """Test conversion of gene-to-gene network to JSON format."""
        g2g_graph = translator.translate_models([readme_model])
        
        # Convert to JSON-serializable format
        g2g_json = {
            "model_info": {
                "id": readme_model.id,
                "title": readme_model.title,
                "taxon": readme_model.taxon,
                "status": readme_model.status
            },
            "nodes": [],
            "edges": []
        }
        
        # Add nodes
        for node, attrs in g2g_graph.nodes(data=True):
            g2g_json["nodes"].append({
                "id": node,
                **attrs
            })
        
        # Add edges
        for source, target, attrs in g2g_graph.edges(data=True):
            g2g_json["edges"].append({
                "source": source,
                "target": target,
                **attrs
            })
        
        # Verify JSON serialization works
        json_str = json.dumps(g2g_json, indent=2)
        assert len(json_str) > 0
        
        # Verify structure
        assert len(g2g_json["nodes"]) == 7
        assert len(g2g_json["edges"]) == 6
        assert g2g_json["model_info"]["id"] == "gomodel:568b0f9600000284"
    
    def test_causal_pathway_structure(self, readme_model, translator):
        """Test that the causal pathway structure is preserved."""
        g2g_graph = translator.translate_models([readme_model])
        
        # Expected causal relationships based on README model
        expected_edges = {
            ("WB:WBGene00006575", "WB:WBGene00003822"),  # tir-1 -> nsy-1
            ("WB:WBGene00006599", "WB:WBGene00012019"),  # tpa-1 -> dkf-2
            ("WB:WBGene00012019", "WB:WBGene00004055"),  # dkf-2 -> pmk-1
            ("WB:WBGene00003822", "WB:WBGene00004758"),  # nsy-1 -> sek-1
            ("WB:WBGene00004758", "WB:WBGene00004055"),  # sek-1 -> pmk-1
            ("WB:WBGene00004055", "WB:WBGene00000223"),  # pmk-1 -> atf-7
        }
        
        actual_edges = set(g2g_graph.edges())
        assert actual_edges == expected_edges
    
    def test_multiple_inputs_to_same_gene(self, readme_model, translator):
        """Test handling of multiple causal inputs to the same gene (pmk-1)."""
        g2g_graph = translator.translate_models([readme_model])
        
        # pmk-1 (WB:WBGene00004055) should have two incoming edges
        pmk1_predecessors = list(g2g_graph.predecessors("WB:WBGene00004055"))
        assert len(pmk1_predecessors) == 2
        
        expected_predecessors = {"WB:WBGene00012019", "WB:WBGene00004758"}  # dkf-2, sek-1
        assert set(pmk1_predecessors) == expected_predecessors
    
    def test_no_duplicate_edges(self, readme_model, translator):
        """Test that no duplicate edges are created."""
        g2g_graph = translator.translate_models([readme_model])
        
        # Count unique edges
        edge_pairs = [(source, target) for source, target in g2g_graph.edges()]
        unique_edges = set(edge_pairs)
        
        # Should be same number (no duplicates)
        assert len(edge_pairs) == len(unique_edges)
    
    def test_all_genes_have_go_annotations(self, readme_model, translator):
        """Test that all edges have some GO term annotations."""
        g2g_graph = translator.translate_models([readme_model])
        
        for source, target, attrs in g2g_graph.edges(data=True):
            # Each edge should have at least some GO terms
            go_attrs = [k for k in attrs.keys() 
                       if k.startswith(("source_gene_", "target_gene_")) 
                       and k.endswith(("_molecular_function", "_biological_process", "_occurs_in"))]
            
            # Should have at least molecular function for both source and target
            assert "source_gene_molecular_function" in attrs
            assert "target_gene_molecular_function" in attrs
            
            # All should have biological process (part_of in this example)
            assert "source_gene_biological_process" in attrs
            assert "target_gene_biological_process" in attrs
    
    def test_readme_example_statistics(self, readme_model, translator):
        """Test that final statistics match expected values."""
        g2g_graph = translator.translate_models([readme_model])
        
        # Original model: 7 activities
        assert len(readme_model.activities) == 7
        
        # Gene-to-gene network: 7 genes, 6 causal relationships
        assert g2g_graph.number_of_nodes() == 7
        assert g2g_graph.number_of_edges() == 6
        
        # All nodes should be gene products
        for node, attrs in g2g_graph.nodes(data=True):
            assert attrs["gene_product"] == node
            assert node.startswith("WB:WBGene")
    
    def test_evidence_collections_basic(self, readme_model, translator):
        """Test that evidence collections are properly included in edges."""
        g2g_graph = translator.translate_models([readme_model])
        
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
    
    def test_molecular_function_evidence_collections(self, readme_model, translator):
        """Test molecular function evidence collections are properly extracted."""
        g2g_graph = translator.translate_models([readme_model])
        
        # Find edge from tir-1 to nsy-1 which should have molecular function evidence
        tir1_to_nsy1_attrs = g2g_graph.edges["WB:WBGene00006575", "WB:WBGene00003822"]
        
        # Check source molecular function evidence
        assert "source_gene_molecular_function_has_reference" in tir1_to_nsy1_attrs
        assert "source_gene_molecular_function_assessed_by" in tir1_to_nsy1_attrs
        
        # Verify reference format
        source_refs = tir1_to_nsy1_attrs["source_gene_molecular_function_has_reference"]
        assert isinstance(source_refs, list)
        assert "PMID:15625192" in source_refs
        
        # Verify evidence code format
        source_codes = tir1_to_nsy1_attrs["source_gene_molecular_function_assessed_by"]
        assert isinstance(source_codes, list)
        assert "ECO:0000314" in source_codes
        
        # Check target molecular function evidence
        assert "target_gene_molecular_function_has_reference" in tir1_to_nsy1_attrs
        assert "target_gene_molecular_function_assessed_by" in tir1_to_nsy1_attrs
        
        target_refs = tir1_to_nsy1_attrs["target_gene_molecular_function_has_reference"]
        assert isinstance(target_refs, list)
        assert "PMID:11751572" in target_refs
    
    def test_biological_process_evidence_collections(self, readme_model, translator):
        """Test biological process evidence collections are properly extracted."""
        g2g_graph = translator.translate_models([readme_model])
        
        # Find edge that should have biological process evidence
        tir1_to_nsy1_attrs = g2g_graph.edges["WB:WBGene00006575", "WB:WBGene00003822"]
        
        # Check source biological process evidence
        assert "source_gene_biological_process_has_reference" in tir1_to_nsy1_attrs
        assert "source_gene_biological_process_assessed_by" in tir1_to_nsy1_attrs
        
        source_refs = tir1_to_nsy1_attrs["source_gene_biological_process_has_reference"]
        assert isinstance(source_refs, list)
        assert "PMID:19837372" in source_refs
        
        source_codes = tir1_to_nsy1_attrs["source_gene_biological_process_assessed_by"]
        assert isinstance(source_codes, list)
        assert "ECO:0000315" in source_codes
    
    def test_occurs_in_evidence_collections(self, readme_model, translator):
        """Test cellular component (occurs_in) evidence collections are properly extracted."""
        g2g_graph = translator.translate_models([readme_model])
        
        # Find edge with occurs_in evidence
        tir1_to_nsy1_attrs = g2g_graph.edges["WB:WBGene00006575", "WB:WBGene00003822"]
        
        # Check source occurs_in evidence
        assert "source_gene_occurs_in_has_reference" in tir1_to_nsy1_attrs
        assert "source_gene_occurs_in_assessed_by" in tir1_to_nsy1_attrs
        
        source_refs = tir1_to_nsy1_attrs["source_gene_occurs_in_has_reference"]
        assert isinstance(source_refs, list)
        assert "PMID:15625192" in source_refs
        
        source_codes = tir1_to_nsy1_attrs["source_gene_occurs_in_assessed_by"]
        assert isinstance(source_codes, list)
        assert "ECO:0000314" in source_codes
    
    def test_causal_predicate_evidence_collections(self, readme_model, translator):
        """Test causal predicate evidence collections are properly extracted."""
        g2g_graph = translator.translate_models([readme_model])
        
        # Find edge with causal association evidence
        tir1_to_nsy1_attrs = g2g_graph.edges["WB:WBGene00006575", "WB:WBGene00003822"]
        
        # Check causal predicate evidence
        assert "causal_predicate_has_reference" in tir1_to_nsy1_attrs
        assert "causal_predicate_assessed_by" in tir1_to_nsy1_attrs
        
        causal_refs = tir1_to_nsy1_attrs["causal_predicate_has_reference"]
        assert isinstance(causal_refs, list)
        assert "PMID:15123841" in causal_refs
        
        causal_codes = tir1_to_nsy1_attrs["causal_predicate_assessed_by"]
        assert isinstance(causal_codes, list)
        assert "ECO:0000315" in causal_codes
    
    def test_contributors_evidence_collections(self, readme_model, translator):
        """Test contributor evidence collections are properly extracted."""
        g2g_graph = translator.translate_models([readme_model])
        
        # Find edge from tir-1 to nsy-1 which should have contributor evidence
        tir1_to_nsy1_attrs = g2g_graph.edges["WB:WBGene00006575", "WB:WBGene00003822"]
        
        # Check for specific contributor evidence that should be present
        expected_contributors = [
            "https://orcid.org/0000-0002-1706-4196",
            "https://orcid.org/0000-0002-3013-9906"
        ]
        
        # Check source molecular function contributors
        assert "source_gene_molecular_function_contributors" in tir1_to_nsy1_attrs
        source_mf_contributors = tir1_to_nsy1_attrs["source_gene_molecular_function_contributors"]
        assert isinstance(source_mf_contributors, list)
        assert "https://orcid.org/0000-0002-1706-4196" in source_mf_contributors
        
        # Check source occurs_in contributors
        assert "source_gene_occurs_in_contributors" in tir1_to_nsy1_attrs
        source_occurs_contributors = tir1_to_nsy1_attrs["source_gene_occurs_in_contributors"]
        assert isinstance(source_occurs_contributors, list)
        assert "https://orcid.org/0000-0002-3013-9906" in source_occurs_contributors
        
        # Check source biological process contributors
        assert "source_gene_biological_process_contributors" in tir1_to_nsy1_attrs
        source_bp_contributors = tir1_to_nsy1_attrs["source_gene_biological_process_contributors"]
        assert isinstance(source_bp_contributors, list)
        assert "https://orcid.org/0000-0002-1706-4196" in source_bp_contributors
        
        # Check causal predicate contributors
        assert "causal_predicate_contributors" in tir1_to_nsy1_attrs
        causal_contributors = tir1_to_nsy1_attrs["causal_predicate_contributors"]
        assert isinstance(causal_contributors, list)
        assert "https://orcid.org/0000-0002-3013-9906" in causal_contributors
    
    def test_evidence_collections_are_lists(self, readme_model, translator):
        """Test that all evidence collection attributes are lists (multivalued)."""
        g2g_graph = translator.translate_models([readme_model])
        
        for source, target, attrs in g2g_graph.edges(data=True):
            for key, value in attrs.items():
                if any(suffix in key for suffix in ['_has_reference', '_assessed_by', '_contributors']):
                    assert isinstance(value, list), f"Evidence collection {key} should be a list, got {type(value)}"
                    assert len(value) > 0, f"Evidence collection {key} should not be empty"
    
    def test_json_serialization_with_evidence_collections(self, readme_model, translator):
        """Test that evidence collections are properly serialized in JSON output."""
        json_output = translator.translate_models_to_json([readme_model])
        
        import json
        result = json.loads(json_output)
        
        # Find edges with evidence collections
        edges_with_evidence = []
        for edge in result["edges"]:
            evidence_attrs = [key for key in edge.keys() 
                             if any(suffix in key for suffix in ['_has_reference', '_assessed_by', '_contributors'])]
            if evidence_attrs:
                edges_with_evidence.append((edge, evidence_attrs))
        
        assert len(edges_with_evidence) > 0, "Should have edges with evidence collections in JSON output"
        
        # Verify structure of evidence collections in JSON
        for edge, evidence_attrs in edges_with_evidence:
            for attr in evidence_attrs:
                assert isinstance(edge[attr], list), f"JSON evidence collection {attr} should be a list"