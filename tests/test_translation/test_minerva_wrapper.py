import json

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


def test_multivalued_input_and_output():
    """Test that activities with multiple inputs and outputs are correctly translated."""
    mw = MinervaWrapper()
    with open(INPUT_DIR / "minerva-633b013300000306.json", "r") as f:
        minerva_object = json.load(f)
    model = mw.minerva_object_to_model(minerva_object)

    cs_activity = next(a for a in model.activities if a.molecular_function.term == "GO:0004108")
    assert len(cs_activity.has_input) == 3
    assert len(cs_activity.has_output) == 2
