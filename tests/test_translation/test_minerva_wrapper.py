import json

import pytest
import yaml

from gocam.translation.minerva_wrapper import MinervaWrapper
from tests import OUTPUT_DIR, INPUT_DIR


# mark as integration test
@pytest.mark.integration
@pytest.mark.parametrize("model_local_id", ["663d668500002178"])
def test_api(model_local_id):
    mw = MinervaWrapper()
    d = mw.internal_object_by_id(model_local_id)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with(open(OUTPUT_DIR / "t.json", "w")) as f:
        json.dump(d, f)
    model = mw.object_by_id(model_local_id)
    assert model is not None


@pytest.mark.parametrize("base_name", ["minerva-example"])
def test_object(base_name):
    mw = MinervaWrapper()
    d = json.load(open(INPUT_DIR / f"{base_name}.json", "r"))
    model = mw.object_from_dict(d)
    assert model is not None
    md = model.model_dump(exclude_unset=True)
    print(yaml.dump(md, sort_keys=False))

