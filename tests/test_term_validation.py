"""
Tests for term validation using linkml-term-validator.

These tests verify that the schema bindings and dynamic enums
correctly validate GO-CAM data against ontology constraints.
"""

from pathlib import Path

import pytest
from linkml.validator import Validator
from linkml.validator.loaders import default_loader_for_file
from linkml.validator.plugins import ValidationPlugin
from linkml_term_validator.plugins import BindingValidationPlugin, DynamicEnumPlugin

SCHEMA_PATH = Path(__file__).parent.parent / "src/gocam/schema/gocam.yaml"


@pytest.fixture
def isolated_cache(tmp_path):
    """Provide an isolated cache directory for each test."""
    return tmp_path / "term_cache"


@pytest.fixture
def binding_validator(isolated_cache):
    """Create a validator with binding validation only."""
    plugins: list[ValidationPlugin] = [
        BindingValidationPlugin(
            oak_adapter_string="sqlite:obo:",
            validate_labels=False,
            strict=True,
            cache_dir=isolated_cache,
        )
    ]
    return Validator(schema=str(SCHEMA_PATH), validation_plugins=plugins)


@pytest.fixture
def full_validator(isolated_cache):
    """Create a validator with both binding and dynamic enum validation."""
    plugins: list[ValidationPlugin] = [
        DynamicEnumPlugin(
            oak_adapter_string="sqlite:obo:",
            cache_dir=isolated_cache,
        ),
        BindingValidationPlugin(
            oak_adapter_string="sqlite:obo:",
            validate_labels=False,
            strict=True,
            cache_dir=isolated_cache,
        ),
    ]
    return Validator(schema=str(SCHEMA_PATH), validation_plugins=plugins)


def make_minimal_model(
    mf_term="GO:0003674",
    bp_term="GO:0008150",
    cc_term="GO:0110165",
    predicate="RO:0002629",
):
    """
    Create a minimal valid GO-CAM model for testing.

    The bindings use `binds_value_of: id` which requires nested objects
    with an `id` field, not plain strings.

    Args:
        mf_term: Molecular function term ID
        bp_term: Biological process term ID
        cc_term: Cellular component term ID
        predicate: Causal predicate ID
    """
    return {
        "id": "gomodel:test-model-001",
        "title": "Test model for validation",
        "taxon": "NCBITaxon:9606",
        "status": "development",
        "activities": [
            {
                "id": "gomodel:test-model-001/activity1",
                "enabled_by": {
                    "type": "EnabledByGeneProductAssociation",
                    "term": "UniProtKB:P12345",
                    "evidence": [],
                    "provenances": [],
                },
                "molecular_function": {
                    "type": "MolecularFunctionAssociation",
                    "term": {"id": mf_term},
                },
                "part_of": {
                    "type": "BiologicalProcessAssociation",
                    "term": {"id": bp_term},
                    "evidence": [],
                    "provenances": [],
                },
                "occurs_in": {
                    "type": "CellularAnatomicalEntityAssociation",
                    "term": {"id": cc_term},
                    "evidence": [],
                    "provenances": [],
                },
                "causal_associations": [
                    {
                        "type": "CausalAssociation",
                        "predicate": {"id": predicate},
                        "downstream_activity": "gomodel:test-model-001/activity2",
                        "evidence": [],
                        "provenances": [],
                    }
                ],
            },
            {
                "id": "gomodel:test-model-001/activity2",
                "enabled_by": {
                    "type": "EnabledByGeneProductAssociation",
                    "term": "UniProtKB:P67890",
                    "evidence": [],
                    "provenances": [],
                },
                "molecular_function": {
                    "type": "MolecularFunctionAssociation",
                    "term": {"id": "GO:0003674"},
                },
            },
        ],
    }


class TestValidModel:
    """Tests with valid GO-CAM models."""

    def test_valid_model_passes_binding_validation(self, binding_validator, tmp_path):
        """A valid model should pass binding validation."""
        import yaml

        model_path = tmp_path / "valid_model.yaml"
        model_path.write_text(yaml.dump(make_minimal_model()))

        loader = default_loader_for_file(model_path)
        report = binding_validator.validate_source(loader, target_class="Model")

        assert len(report.results) == 0, f"Expected no errors, got: {report.results}"


class TestInvalidPredicates:
    """Tests for invalid causal predicates."""

    def test_nonexistent_predicate_fails(self, binding_validator, tmp_path):
        """A predicate ID that doesn't exist should fail validation."""
        import yaml

        model = make_minimal_model(predicate="RO:9999999")
        model_path = tmp_path / "bad_predicate.yaml"
        model_path.write_text(yaml.dump(model))

        loader = default_loader_for_file(model_path)
        report = binding_validator.validate_source(loader, target_class="Model")

        assert len(report.results) > 0, (
            "Expected validation errors for nonexistent predicate"
        )
        error_messages = [r.message for r in report.results]
        assert any("RO:9999999" in msg for msg in error_messages), (
            f"Expected error mentioning RO:9999999, got: {error_messages}"
        )

    def test_wrong_predicate_type_fails(self, binding_validator, tmp_path):
        """Using a non-causal RO relation should fail binding validation."""
        import yaml

        # RO:0000052 is 'inheres in' - not a causal predicate
        model = make_minimal_model(predicate="RO:0000052")
        model_path = tmp_path / "wrong_predicate.yaml"
        model_path.write_text(yaml.dump(model))

        loader = default_loader_for_file(model_path)
        report = binding_validator.validate_source(loader, target_class="Model")

        assert len(report.results) > 0, (
            "Expected validation errors for wrong predicate type"
        )


class TestInvalidOntologyTerms:
    """Tests for invalid ontology terms in associations."""

    @pytest.mark.integration
    def test_nonexistent_go_term_fails(self, full_validator, tmp_path):
        """A GO term ID that doesn't exist should fail validation."""
        import yaml

        model = make_minimal_model(mf_term="GO:9999999")
        model_path = tmp_path / "bad_go_term.yaml"
        model_path.write_text(yaml.dump(model))

        loader = default_loader_for_file(model_path)
        report = full_validator.validate_source(loader, target_class="Model")

        assert len(report.results) > 0, (
            "Expected validation errors for nonexistent GO term"
        )
        error_messages = [r.message for r in report.results]
        assert any("GO:9999999" in msg for msg in error_messages), (
            f"Expected error mentioning GO:9999999, got: {error_messages}"
        )

    @pytest.mark.integration
    def test_cc_term_in_bp_slot_fails(self, full_validator, tmp_path):
        """Using a cellular component term where biological process is expected should fail."""
        import yaml

        # GO:0005634 is 'nucleus' - a cellular component, not a biological process
        model = make_minimal_model(bp_term="GO:0005634")
        model_path = tmp_path / "cc_in_bp_slot.yaml"
        model_path.write_text(yaml.dump(model))

        loader = default_loader_for_file(model_path)
        report = full_validator.validate_source(loader, target_class="Model")

        assert len(report.results) > 0, (
            "Expected validation errors for CC term in BP slot"
        )
        error_messages = [r.message for r in report.results]
        assert any("GO:0005634" in msg for msg in error_messages), (
            f"Expected error mentioning GO:0005634, got: {error_messages}"
        )

    @pytest.mark.integration
    def test_bp_term_in_cc_slot_fails(self, full_validator, tmp_path):
        """Using a biological process term where cellular component is expected should fail."""
        import yaml

        # GO:0006915 is 'apoptotic process' - a BP, not a CC
        model = make_minimal_model(cc_term="GO:0006915")
        model_path = tmp_path / "bp_in_cc_slot.yaml"
        model_path.write_text(yaml.dump(model))

        loader = default_loader_for_file(model_path)
        report = full_validator.validate_source(loader, target_class="Model")

        assert len(report.results) > 0, (
            "Expected validation errors for BP term in CC slot"
        )

    @pytest.mark.integration
    def test_mf_term_in_bp_slot_fails(self, full_validator, tmp_path):
        """Using a molecular function term where biological process is expected should fail."""
        import yaml

        # GO:0003824 is 'catalytic activity' - an MF, not a BP
        model = make_minimal_model(bp_term="GO:0003824")
        model_path = tmp_path / "mf_in_bp_slot.yaml"
        model_path.write_text(yaml.dump(model))

        loader = default_loader_for_file(model_path)
        report = full_validator.validate_source(loader, target_class="Model")

        assert len(report.results) > 0, (
            "Expected validation errors for MF term in BP slot"
        )


class TestExistingDataValidation:
    """Tests that validate existing data files."""

    @pytest.mark.integration
    def test_existing_model_passes(self, full_validator):
        """An existing model from the test fixtures should pass validation."""
        model_path = Path(__file__).parent / "input/Model-63f809ec00000701.yaml"
        if not model_path.exists():
            pytest.skip("Test model file not found")

        loader = default_loader_for_file(model_path)
        report = full_validator.validate_source(loader, target_class="Model")

        assert len(report.results) == 0, (
            f"Expected no errors for existing model, got: {report.results}"
        )
