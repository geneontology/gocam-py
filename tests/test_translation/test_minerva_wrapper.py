import json

import pytest

from gocam.translation.minerva_wrapper import MinervaWrapper
from tests import INPUT_DIR

ENABLED_BY = "RO:0002333"


# mark as integration test
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
