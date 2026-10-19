"""Data test."""

import glob
import os

import pytest
from linkml.validator import Validator
from linkml.validator.plugins import JsonschemaValidationPlugin
from linkml_runtime.loaders import yaml_loader

from gocam.datamodel import MAIN_SCHEMA_PATH, Model
from tests import EXAMPLES_DIR

EXAMPLE_FILES = glob.glob(os.path.join(EXAMPLES_DIR, "*.yaml"))


@pytest.fixture(scope="module")
def validator():
    return Validator(
        schema=MAIN_SCHEMA_PATH,
        validation_plugins=[JsonschemaValidationPlugin(closed=True)],
    )


# linkml-linkml_runtime.loaders.yaml_loader still uses load_obj. Disable capturing this warning.
@pytest.mark.filterwarnings("ignore::pydantic.warnings.PydanticDeprecatedSince20")
def test_data():
    """Data test."""
    for path in EXAMPLE_FILES:
        obj = yaml_loader.load(path, target_class=Model)
        assert obj


@pytest.mark.parametrize(
    "uses_placeholder_term,has_part_present,valid",
    [
        (True, True, True),  # Placeholder term with has_part present is valid
        (True, False, False),  # Placeholder term with has_part missing is invalid
        (False, True, True),  # Non-placeholder term with has_part present is valid
        (False, False, True),  # Non-placeholder term with has_part missing is valid
    ],
)
def test_protein_complex_association_rules(
    uses_placeholder_term, has_part_present, valid, validator
):
    """Test the rules on the EnabledByProteinComplexAssociation class."""

    instance = {
        "type": "EnabledByProteinComplexAssociation",
        "term": "GO:0032991" if uses_placeholder_term else "GO:0097141",
    }
    if has_part_present:
        instance["has_part"] = [
            {"type": "ProteinComplexHasPartAssociation", "term": "UniProtKB:O75425"},
            {"type": "ProteinComplexHasPartAssociation", "term": "UniProtKB:Q9UJG1"},
        ]
    report = validator.validate(
        instance, target_class="EnabledByProteinComplexAssociation"
    )
    if valid:
        assert not report.results
    else:
        assert len(report.results) == 1
        assert "has_part" in report.results[0].message
