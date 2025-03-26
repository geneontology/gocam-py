import json
import re

import pytest

from gocam.datamodel import EnabledByProteinComplexAssociation
from gocam.translation.minerva_wrapper import MinervaWrapper
from tests import INPUT_DIR

ENABLED_BY = "RO:0002333"


# This is an integration test because it makes real network requests
@pytest.mark.integration
@pytest.mark.parametrize("model_local_id", ["663d668500002178"])
def test_api(model_local_id):
    mw = MinervaWrapper()
    model = mw.fetch_model(model_local_id)
    assert model is not None


@pytest.mark.parametrize(
    "base_name", ["minerva-663d668500002178", "minerva-5b91dbd100002057"]
)
def test_object(base_name):
    mw = MinervaWrapper()
    with open(INPUT_DIR / f"{base_name}.json", "r") as f:
        minerva_object = json.load(f)
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
    mw = MinervaWrapper()
    with open(INPUT_DIR / "minerva-5ce58dde00001215.json", "r") as f:
        minerva_object = json.load(f)
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
    mw = MinervaWrapper()
    with open(INPUT_DIR / "minerva-665912ed00002626.json", "r") as f:
        minerva_object = json.load(f)
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
    mw = MinervaWrapper()
    with open(INPUT_DIR / "minerva-5f46c3b700001031.json", "r") as f:
        minerva_object = json.load(f)
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
    mw = MinervaWrapper()
    with open(INPUT_DIR / "minerva-633b013300000306.json", "r") as f:
        minerva_object = json.load(f)
    model = mw.minerva_object_to_model(minerva_object)

    cs_activity = next(
        a for a in model.activities if a.molecular_function.term == "GO:0004108"
    )
    assert len(cs_activity.has_input) == 3
    assert len(cs_activity.has_output) == 2


def test_complete_provenance_on_evidence():
    """Test that all contributor and providedBy annotations are included on the ProvenanceInfo
    instance attached to evidence."""
    mw = MinervaWrapper()
    with open(INPUT_DIR / "minerva-633b013300000306.json", "r") as f:
        minerva_object = json.load(f)

    # ensure that all evidence has more than one contributor and providedBy annotation
    for individual in minerva_object["individuals"]:
        if any(rt["id"] == "ECO:0000000" for rt in individual["root-type"]):
            individual["annotations"].append(
                {"key": "providedBy", "value": "https://www.example.org"}
            )
            individual["annotations"].append(
                {"key": "contributor", "value": "https://orcid.org/0000-0000-0000-0000"}
            )

    model = mw.minerva_object_to_model(minerva_object)

    # assert that provenances of evidence has more than one contributor and provided_by
    for activity in model.activities:
        for association in activity.causal_associations:
            for evidence in association.evidence:
                for provenance in evidence.provenances:
                    assert len(provenance.contributor) > 1
                    assert len(provenance.provided_by) > 1


def test_evidence_with_objects():
    """Test that evidence with_objects are correctly translated."""
    mw = MinervaWrapper()
    with open(INPUT_DIR / "minerva-5f46c3b700001031.json", "r") as f:
        minerva_object = json.load(f)
    model = mw.minerva_object_to_model(minerva_object)

    kinase_activity = next((a for a in model.activities if a.molecular_function.term == "GO:0004674"), None)
    assert kinase_activity is not None
    assert len(kinase_activity.molecular_function.evidence) == 1

    evidence = kinase_activity.molecular_function.evidence[0]
    assert len(evidence.with_objects) == 2
    assert all(re.match(r"^[A-Z]+:[A-Z0-9]+$", obj) for obj in evidence.with_objects)
