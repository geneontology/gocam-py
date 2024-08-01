"""Data test."""

import glob
import os

import pytest
from linkml_runtime.loaders import yaml_loader

from gocam.datamodel import Model
from tests import EXAMPLES_DIR

EXAMPLE_FILES = glob.glob(os.path.join(EXAMPLES_DIR, "*.yaml"))


# linkml-linkml_runtime.loaders.yaml_loader still uses load_obj. Disable capturing this warning.
@pytest.mark.filterwarnings("ignore::pydantic.warnings.PydanticDeprecatedSince20")
def test_data():
    """Data test."""
    for path in EXAMPLE_FILES:
        obj = yaml_loader.load(path, target_class=Model)
        assert obj
