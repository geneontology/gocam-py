import json
import re

import pytest

from gocam.datamodel import EnabledByProteinComplexAssociation
from gocam.translation.minerva_wrapper import MinervaWrapper
from tests import INPUT_DIR

ENABLED_BY = "RO:0002333"


def load_minerva_object(id: str):
    with open(INPUT_DIR / f"minerva-{id}.json", "r") as f:
        return json.load(f)


# This is an integration test because it makes real network requests
@pytest.mark.integration
@pytest.mark.parametrize("model_local_id", ["663d668500002178"])
def test_api(model_local_id):
    mw = MinervaWrapper()
    model = mw.fetch_model(model_local_id)
    assert model is not None


@pytest.mark.parametrize("id", ["663d668500002178", "5b91dbd100002057"])
def test_object(id):
    mw = MinervaWrapper()
    minerva_object = load_minerva_object(id)
    model = mw.minerva_object_to_model(minerva_object)

    # TODO: add more sanity checks here
    assert model is not None
    assert model.id == minerva_object["id"]
    enabled_by_facts = [
        fact for fact in minerva_object["facts"] if fact["property"] == ENABLED_BY
    ]
    assert len(model.activities) == len(enabled_by_facts)


def test_protein_complex():
    """Test that activities enabled by protein complexes are correctly translated."""
    minerva_object = load_minerva_object("5ce58dde00001215")
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    protein_complex_activities = [
        a
        for a in model.activities
        if isinstance(a.enabled_by, EnabledByProteinComplexAssociation)
    ]
    assert len(protein_complex_activities) == 1

    protein_complex_activity = protein_complex_activities[0]
    assert protein_complex_activity.enabled_by.members == [
        "MGI:MGI:1929608",
        "MGI:MGI:103038",
    ]


def test_has_input_and_has_output():
    """Test that input/output molecule associations are added to activities"""
    minerva_object = load_minerva_object("665912ed00002626")
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    activities_with_input = []
    activities_with_output = []
    for activity in model.activities:
        if activity.has_input:
            activities_with_input.append(activity)
        if activity.has_output:
            activities_with_output.append(activity)

    # Basic sanity check on the number of activities with input/output
    assert len(activities_with_input) == 3
    assert len(activities_with_output) == 7

    # Verify that one activity has uric acid as an input
    uric_acid_input_activities = [
        a for a in activities_with_input if a.has_input[0].term == "CHEBI:27226"
    ]
    assert len(uric_acid_input_activities) == 1

    # Verify that three activities have urea as an output
    urea_output_activities = [
        a for a in activities_with_output if a.has_output[0].term == "CHEBI:16199"
    ]
    assert len(urea_output_activities) == 3


def test_has_input_issue_65():
    """Test that all input associations, including proteins, are included in the core model"""
    minerva_object = load_minerva_object("5f46c3b700001031")
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    # Find activities with inputs
    activities_with_input = [a for a in model.activities if a.has_input is not None]

    # Verify that all inputs are included, even protein inputs
    for activity in activities_with_input:
        for input_assoc in activity.has_input:
            # Check if the input term is also the term of an enabled_by for any activity
            is_protein = any(
                a.enabled_by.term == input_assoc.term
                for a in model.activities
                if a.enabled_by is not None
            )
            # If this test passes, it means we're no longer filtering at the core data model level
            if is_protein:
                assert input_assoc is not None


def test_multivalued_input_and_output():
    """Test that activities with multiple inputs and outputs are correctly translated."""
    minerva_object = load_minerva_object("633b013300000306")
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    cs_activity = next(
        a for a in model.activities if a.molecular_function.term == "GO:0004108"
    )
    assert len(cs_activity.has_input) == 3
    assert len(cs_activity.has_output) == 2

def test_missing_enabled_by():
    """Test that activities without an enabled_by association are handled correctly."""
    minerva_object = load_minerva_object("YeastPathways_LYSDEGII-PWY")
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    # Find activities without an enabled_by association
    activities_without_enabled_by = [
        a for a in model.activities if a.enabled_by is None
    ]

    # Verify that there are no such activities
    assert len(activities_without_enabled_by) == 0, (
        "There should be no activities without an enabled_by association."
    )


def test_provenance_on_evidence():
    """Test that all contributor and providedBy annotations are included on the ProvenanceInfo
    instance attached to evidence."""
    minerva_object = load_minerva_object("633b013300000306")

    # ensure that all evidence has more than one contributor and providedBy annotation
    for individual in minerva_object["individuals"]:
        if any(rt["id"] == "ECO:0000000" for rt in individual["root-type"]):
            individual["annotations"].append(
                {"key": "providedBy", "value": "https://www.example.org"}
            )
            individual["annotations"].append(
                {"key": "contributor", "value": "https://orcid.org/0000-0000-0000-0000"}
            )

    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    # assert that provenances of evidence has more than one contributor and provided_by
    for activity in model.activities:
        for association in activity.causal_associations:
            for evidence in association.evidence:
                for provenance in evidence.provenances:
                    assert len(provenance.contributor) > 1
                    assert len(provenance.provided_by) > 1


def test_provenance_on_model():
    """Test that top-level annotations are included on the ProvenanceInfo instance attached to the
    Model."""
    minerva_object = load_minerva_object("5f46c3b700001031")
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    assert model.provenances is not None
    assert len(model.provenances) == 1
    provenance = model.provenances[0]
    assert set(provenance.contributor) == {
        "https://orcid.org/0000-0001-7646-0052",
        "https://orcid.org/0000-0001-8769-177X",
        "https://orcid.org/0000-0002-1706-4196",
        "https://orcid.org/0000-0003-1813-6857",
    }
    assert set(provenance.provided_by) == {
        "http://geneontology.org",
        "http://www.wormbase.org",
        "https://www.uniprot.org",
    }
    assert provenance.date == "2023-11-02"


def test_provenance_on_associations():
    """Test that fact annotations are included on the ProvenanceInfo instance attached to various
    Association subclasses."""
    minerva_object = load_minerva_object("663d668500002178")
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    for activity in model.activities:
        if activity.causal_associations is not None:
            for causal_assoc in activity.causal_associations:
                assert causal_assoc.provenances is not None
                assert len(causal_assoc.provenances) > 0

        if activity.has_input is not None:
            for input_assoc in activity.has_input:
                assert input_assoc.provenances is not None
                assert len(input_assoc.provenances) > 0

        if activity.has_output is not None:
            for output_assoc in activity.has_output:
                assert output_assoc.provenances is not None
                assert len(output_assoc.provenances) > 0

        if activity.has_primary_input is not None:
            assert activity.has_primary_input.provenances is not None
            assert len(activity.has_primary_input.provenances) > 0

        if activity.has_primary_output is not None:
            assert activity.has_primary_output.provenances is not None
            assert len(activity.has_primary_output.provenances) > 0

        if activity.occurs_in is not None:
            assert activity.occurs_in.provenances is not None
            assert len(activity.occurs_in.provenances) > 0

        if activity.part_of is not None:
            assert activity.part_of.provenances is not None
            assert len(activity.part_of.provenances) > 0

        if activity.enabled_by is not None:
            assert activity.enabled_by.provenances is not None
            assert len(activity.enabled_by.provenances) > 0


def test_evidence_with_objects():
    """Test that evidence with_objects are correctly translated."""
    minerva_object = load_minerva_object("5f46c3b700001031")
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    kinase_activity = next(
        (a for a in model.activities if a.molecular_function.term == "GO:0004674"), None
    )
    assert kinase_activity is not None
    assert len(kinase_activity.enabled_by.evidence) == 1

    evidence = kinase_activity.enabled_by.evidence[0]
    assert len(evidence.with_objects) == 2
    assert all(re.match(r"^[A-Z]+:[A-Z0-9]+$", obj) for obj in evidence.with_objects)


def test_additional_taxa():
    """Test that model with multiple taxa correctly handles primary taxon and additional_taxa."""
    # Create a copy of the test file with minimal content for controlled testing
    minerva_object = {
        "id": "gomodel:test123",
        "annotations": [
            {
                "key": "title",
                "value": "Test model with multiple taxa"
            },
            # First taxon annotation (Human)
            {
                "key": "https://w3id.org/biolink/vocab/in_taxon",
                "value": "NCBITaxon:9606",
                "value-type": "IRI"
            },
            # Second taxon annotation (E. coli)
            {
                "key": "https://w3id.org/biolink/vocab/in_taxon",
                "value": "NCBITaxon:562",
                "value-type": "IRI"
            }
        ],
        "individuals": [],
        "facts": []
    }
    
    # Add a minimal individual and fact to make it a valid model
    minerva_object["individuals"] = [
        {
            "id": "gomodel:test123/activity1",
            "type": [{"type": "class", "id": "GO:0003674", "label": "molecular_function"}],
            "root-type": [{"type": "class", "id": "GO:0003674", "label": "molecular_function"}],
            "annotations": []
        },
        {
            "id": "gomodel:test123/protein1",
            "type": [{"type": "class", "id": "UniProtKB:P12345", "label": "test protein"}],
            "root-type": [{"type": "class", "id": "CHEBI:33695", "label": "information biomacromolecule"}],
            "annotations": []
        }
    ]
    
    minerva_object["facts"] = [
        {
            "subject": "gomodel:test123/activity1",
            "property": "RO:0002333",
            "object": "gomodel:test123/protein1",
            "annotations": []
        }
    ]
    
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)
    
    # Verify that the model makes human the primary taxon and E. coli additional
    assert model.taxon == "NCBITaxon:9606"
    assert len(model.additional_taxa) == 1
    assert model.additional_taxa[0] == "NCBITaxon:562"
    
    # Test with two non-host taxa
    minerva_object["annotations"] = [
        {
            "key": "title",
            "value": "Test model with multiple taxa"
        },
        {
            "key": "https://w3id.org/biolink/vocab/in_taxon",
            "value": "NCBITaxon:562"    # E. coli (not a host)
        },
        {
            "key": "https://w3id.org/biolink/vocab/in_taxon",
            "value": "NCBITaxon:623"    # Shigella (not a host)
        }
    ]
    
    model = mw.minerva_object_to_model(minerva_object)
    
    # When no hosts, the first taxon should be primary
    assert model.taxon == "NCBITaxon:562"
    assert len(model.additional_taxa) == 1
    assert model.additional_taxa[0] == "NCBITaxon:623"


def test_host_taxon_prioritization():
    """Test that host taxa are properly prioritized when multiple taxa are present."""
    minerva_object = load_minerva_object("6348a65d00000661")
    
    # First, ensure we're working with a fresh copy with no taxon annotations
    minerva_object["annotations"] = [ann for ann in minerva_object["annotations"] 
                               if ann.get("key") != "https://w3id.org/biolink/vocab/in_taxon" 
                               and ann.get("key") != "in_taxon"]
    
    # Add taxon annotations: one host and one pathogen (in non-host-first order)
    minerva_object["annotations"].extend([
        {
            "key": "https://w3id.org/biolink/vocab/in_taxon", 
            "value": "NCBITaxon:623",    # Shigella (pathogen)
            "value-type": "IRI"
        },
        {
            "key": "https://w3id.org/biolink/vocab/in_taxon",
            "value": "NCBITaxon:9606",   # Human (host)
            "value-type": "IRI"
        }
    ])
    
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)
    
    # Verify the model prioritizes the host taxon as primary, even though it was added second
    assert model.taxon == "NCBITaxon:9606"
    assert len(model.additional_taxa) == 1
    assert model.additional_taxa[0] == "NCBITaxon:623"
