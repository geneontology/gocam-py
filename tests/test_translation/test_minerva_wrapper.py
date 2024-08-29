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
