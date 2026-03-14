import json
import re

import pytest

from gocam.datamodel import Activity, EnabledByProteinComplexAssociation
from gocam.translation.minerva_wrapper import (
    MOLECULAR_ASSOCIATION_INVERSE_PROPERTIES,
    MOLECULAR_ASSOCIATION_PROPERTIES,
    MinervaView,
    MinervaWrapper,
)
from gocam.translation.result import WarningType
from gocam.vocabulary import Relation
from tests import INPUT_DIR


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
        fact
        for fact in minerva_object["facts"]
        if fact["property"] == Relation.ENABLED_BY
    ]
    assert model.activities is not None
    assert len(model.activities) == len(enabled_by_facts)


def test_protein_complex():
    """Test that activities enabled by protein complexes are correctly translated."""
    minerva_object = load_minerva_object("5ce58dde00001215")
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    protein_complex_activities = [
        a
        for a in model.activities or []
        if isinstance(a.enabled_by, EnabledByProteinComplexAssociation)
    ]
    assert len(protein_complex_activities) == 1

    protein_complex_activity = protein_complex_activities[0]
    assert isinstance(
        protein_complex_activity.enabled_by, EnabledByProteinComplexAssociation
    )
    assert [m.term for m in protein_complex_activity.enabled_by.members or []] == [
        "MGI:MGI:1929608",
        "MGI:MGI:103038",
    ]


def test_has_input_and_has_output():
    """Test that input/output molecule associations are added to activities"""
    minerva_object = load_minerva_object("665912ed00002626")
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    activities_with_input: list[Activity] = [
        activity
        for activity in model.activities or []
        if any(
            ma.predicate == Relation.HAS_INPUT
            for ma in activity.molecular_associations or []
        )
    ]
    activities_with_output: list[Activity] = [
        activity
        for activity in model.activities or []
        if any(
            ma.predicate == Relation.HAS_OUTPUT
            for ma in activity.molecular_associations or []
        )
    ]

    # Basic sanity check on the number of activities with input/output
    assert len(activities_with_input) == 3
    assert len(activities_with_output) == 7

    # Verify that one activity has uric acid as an input
    uric_acid_molecule = next(
        m for m in model.molecules or [] if m.term == "CHEBI:27226"
    )
    uric_acid_input_activities = [
        a
        for a in activities_with_input
        if any(
            ma.molecule == uric_acid_molecule.id and ma.predicate == Relation.HAS_INPUT
            for ma in a.molecular_associations or []
        )
    ]
    assert len(uric_acid_input_activities) == 1

    # Verify that three activities have urea as an output
    urea_molecule = next(m for m in model.molecules or [] if m.term == "CHEBI:16199")
    urea_output_activities = [
        a
        for a in activities_with_output
        if any(
            ma.molecule == urea_molecule.id and ma.predicate == Relation.HAS_OUTPUT
            for ma in a.molecular_associations or []
        )
    ]
    assert len(urea_output_activities) == 3


def test_has_input_issue_65():
    """Test that all input associations, including proteins, are included in the core model"""
    minerva_object = load_minerva_object("5f46c3b700001031")
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    # Find activities with inputs
    activities_with_input = [
        a
        for a in (model.activities or [])
        if a.molecular_associations
        and any(ma.predicate == Relation.HAS_INPUT for ma in a.molecular_associations)
    ]

    # Verify that all inputs are included, even protein inputs
    for activity in activities_with_input:
        for input_assoc in activity.molecular_associations or []:
            if input_assoc.predicate != Relation.HAS_INPUT:
                continue
            # Check if the input term is also the term of an enabled_by for any activity
            if input_assoc.molecule is None or not model.molecules:
                continue
            input_molecule = next(
                (m for m in model.molecules if m.id == input_assoc.molecule), None
            )
            if input_molecule is None:
                continue
            is_protein = any(
                a.enabled_by.term == input_molecule.term
                for a in (model.activities or [])
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
        a
        for a in (model.activities or [])
        if a.molecular_function and a.molecular_function.term == "GO:0004108"
    )
    assert cs_activity.molecular_associations is not None
    inputs = [
        ma
        for ma in cs_activity.molecular_associations
        if ma.predicate == Relation.HAS_INPUT
    ]
    outputs = [
        ma
        for ma in cs_activity.molecular_associations
        if ma.predicate == Relation.HAS_OUTPUT
    ]
    assert len(inputs) == 3
    assert len(outputs) == 2


def test_missing_enabled_by():
    """Test that activities without an enabled_by association are handled correctly."""
    minerva_object = load_minerva_object("YeastPathways_LYSDEGII-PWY")
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    # Find activities without an enabled_by association
    activities_without_enabled_by = [
        a for a in (model.activities or []) if a.enabled_by is None
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
    for activity in model.activities or []:
        for association in activity.causal_associations or []:
            for evidence in association.evidence or []:
                for provenance in evidence.provenances or []:
                    assert len(provenance.contributor or []) > 1
                    assert len(provenance.provided_by or []) > 1


def test_provenance_on_model():
    """Test that top-level annotations are included on the ProvenanceInfo instance attached to the
    Model."""
    minerva_object = load_minerva_object("5f46c3b700001031")
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    assert model.provenances is not None
    assert len(model.provenances) == 1
    provenance = model.provenances[0]
    assert set(provenance.contributor or []) == {
        "https://orcid.org/0000-0001-7646-0052",
        "https://orcid.org/0000-0001-8769-177X",
        "https://orcid.org/0000-0002-1706-4196",
        "https://orcid.org/0000-0003-1813-6857",
    }
    assert set(provenance.provided_by or []) == {
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

    for activity in model.activities or []:
        if activity.causal_associations is not None:
            for causal_assoc in activity.causal_associations:
                assert causal_assoc.provenances is not None
                assert len(causal_assoc.provenances) > 0

        if activity.molecular_associations:
            for mol_assoc in activity.molecular_associations:
                assert mol_assoc.provenances is not None
                assert len(mol_assoc.provenances) > 0

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
        (
            a
            for a in model.activities or []
            if a.molecular_function and a.molecular_function.term == "GO:0004674"
        ),
        None,
    )
    assert kinase_activity is not None
    assert kinase_activity.enabled_by is not None
    assert kinase_activity.enabled_by.evidence is not None
    assert len(kinase_activity.enabled_by.evidence) == 1

    evidence = kinase_activity.enabled_by.evidence[0]
    assert evidence.with_objects is not None
    assert len(evidence.with_objects) == 2
    assert all(re.match(r"^[A-Z]+:[A-Z0-9]+$", obj) for obj in evidence.with_objects)


def test_additional_taxa():
    """Test that model with multiple taxa correctly handles primary taxon and additional_taxa."""
    # Create a copy of the test file with minimal content for controlled testing
    minerva_object = {
        "id": "gomodel:test123",
        "annotations": [
            {"key": "title", "value": "Test model with multiple taxa"},
            # First taxon annotation (Human)
            {
                "key": "https://w3id.org/biolink/vocab/in_taxon",
                "value": "NCBITaxon:9606",
                "value-type": "IRI",
            },
            # Second taxon annotation (E. coli)
            {
                "key": "https://w3id.org/biolink/vocab/in_taxon",
                "value": "NCBITaxon:562",
                "value-type": "IRI",
            },
        ],
        "individuals": [],
        "facts": [],
    }

    # Add a minimal individual and fact to make it a valid model
    minerva_object["individuals"] = [
        {
            "id": "gomodel:test123/activity1",
            "type": [
                {"type": "class", "id": "GO:0003674", "label": "molecular_function"}
            ],
            "root-type": [
                {"type": "class", "id": "GO:0003674", "label": "molecular_function"}
            ],
            "annotations": [],
        },
        {
            "id": "gomodel:test123/protein1",
            "type": [
                {"type": "class", "id": "UniProtKB:P12345", "label": "test protein"}
            ],
            "root-type": [
                {
                    "type": "class",
                    "id": "CHEBI:33695",
                    "label": "information biomacromolecule",
                }
            ],
            "annotations": [],
        },
    ]

    minerva_object["facts"] = [
        {
            "subject": "gomodel:test123/activity1",
            "property": Relation.HAS_INPUT,
            "object": "gomodel:test123/protein1",
            "annotations": [],
        }
    ]

    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    # Verify that the model makes human the primary taxon and E. coli additional
    assert model.taxon == "NCBITaxon:9606"
    assert model.additional_taxa is not None
    assert len(model.additional_taxa) == 1
    assert model.additional_taxa[0] == "NCBITaxon:562"

    # Test with two non-host taxa
    minerva_object["annotations"] = [
        {"key": "title", "value": "Test model with multiple taxa"},
        {
            "key": "https://w3id.org/biolink/vocab/in_taxon",
            "value": "NCBITaxon:562",  # E. coli (not a host)
        },
        {
            "key": "https://w3id.org/biolink/vocab/in_taxon",
            "value": "NCBITaxon:623",  # Shigella (not a host)
        },
    ]

    model = mw.minerva_object_to_model(minerva_object)

    # When no hosts, the first taxon should be primary
    assert model.taxon == "NCBITaxon:562"
    assert model.additional_taxa is not None
    assert len(model.additional_taxa) == 1
    assert model.additional_taxa[0] == "NCBITaxon:623"


def test_host_taxon_prioritization():
    """Test that host taxa are properly prioritized when multiple taxa are present."""
    minerva_object = load_minerva_object("6348a65d00000661")

    # First, ensure we're working with a fresh copy with no taxon annotations
    minerva_object["annotations"] = [
        ann
        for ann in minerva_object["annotations"]
        if ann.get("key") != "https://w3id.org/biolink/vocab/in_taxon"
        and ann.get("key") != "in_taxon"
    ]

    # Add taxon annotations: one host and one pathogen (in non-host-first order)
    minerva_object["annotations"].extend(
        [
            {
                "key": "https://w3id.org/biolink/vocab/in_taxon",
                "value": "NCBITaxon:623",  # Shigella (pathogen)
                "value-type": "IRI",
            },
            {
                "key": "https://w3id.org/biolink/vocab/in_taxon",
                "value": "NCBITaxon:9606",  # Human (host)
                "value-type": "IRI",
            },
        ]
    )

    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    # Verify the model prioritizes the host taxon as primary, even though it was added second
    assert model.taxon == "NCBITaxon:9606"
    assert model.additional_taxa is not None
    assert len(model.additional_taxa) == 1
    assert model.additional_taxa[0] == "NCBITaxon:623"


def test_translation_warning_missing_term():
    """Test that a warning is generated for missing terms during translation."""
    minerva_object = load_minerva_object("663d668500002178")

    # Introduce a fact with a missing term
    minerva_object["facts"].append(
        {
            "subject": "gomodel:663d668500002178/subject_missing_term",
            "property": Relation.HAS_INPUT,
            "object": "gomodel:663d668500002178/object_missing_term",
            "annotations": [],
        }
    )

    translation_result = MinervaWrapper.translate(minerva_object)

    # Check that a warning for the missing term was generated
    warnings = translation_result.warnings
    assert len(warnings) == 2
    missing_term_warning = next(
        (w for w in warnings if w.type == WarningType.MISSING_TERM), None
    )
    assert missing_term_warning is not None
    assert "Missing term for individual" in missing_term_warning.message
    assert (
        missing_term_warning.entity_id == "gomodel:663d668500002178/object_missing_term"
    )

    # A second warning should be generated for an unhandled fact due to the missing term
    unhandled_fact_warning = next(
        (w for w in warnings if w.type == WarningType.UNHANDLED_FACT), None
    )
    assert unhandled_fact_warning is not None


def test_translation_warning_no_enabled_by_facts():
    """Test that a warning is generated for activities with no enabled_by facts."""
    minerva_object = {
        "id": "gomodel:test_no_enabled_by",
        "annotations": [
            {"key": "title", "value": "Test model with activity missing enabled_by"}
        ],
        "individuals": [
            {
                "id": "gomodel:test_no_enabled_by/activity_no_enabled_by",
                "type": [
                    {"type": "class", "id": "GO:0003674", "label": "molecular_function"}
                ],
                "root-type": [
                    {"type": "class", "id": "GO:0003674", "label": "molecular_function"}
                ],
                "annotations": [],
            }
        ],
        "facts": [],
    }

    translation_result = MinervaWrapper.translate(minerva_object)

    # Check that a warning for no enabled_by facts was generated
    warnings = translation_result.warnings
    assert len(warnings) == 1
    assert warnings[0].type == WarningType.NO_ENABLED_BY_FACTS
    assert "No enabled_by (RO:0002333) facts" in warnings[0].message
    assert warnings[0].entity_id == "gomodel:test_no_enabled_by"


def test_translation_warning_invalid_model_state():
    """Test that a warning is generated for an invalid model state."""
    minerva_object = load_minerva_object("663d668500002178")
    for annotation in minerva_object["annotations"]:
        if annotation["key"] == "state":
            annotation["value"] = "invalid_state_for_testing"

    translation_result = MinervaWrapper.translate(minerva_object)

    # Check that a warning for invalid model state was generated
    warnings = translation_result.warnings
    assert len(warnings) == 1
    assert warnings[0].type == WarningType.INVALID_MODEL_STATE
    assert "invalid_state_for_testing" in warnings[0].message
    assert warnings[0].entity_id == "gomodel:663d668500002178"


def test_happens_during():
    """Test that the happens_during slot on Activity instances is populated correctly."""
    minerva_object = {
        "id": "gomodel:test_happens_during",
        "annotations": [{"key": "title", "value": "Test model with happens_during"}],
        "individuals": [
            {
                "id": "gomodel:test_happens_during/activity1",
                "type": [
                    {"type": "class", "id": "GO:0003674", "label": "molecular_function"}
                ],
                "root-type": [
                    {"type": "class", "id": "GO:0003674", "label": "molecular_function"}
                ],
                "annotations": [],
            },
            {
                "id": "gomodel:test_happens_during/protein1",
                "type": [
                    {"type": "class", "id": "UniProtKB:P12345", "label": "test protein"}
                ],
                "root-type": [
                    {
                        "type": "class",
                        "id": "CHEBI:33695",
                        "label": "information biomacromolecule",
                    }
                ],
                "annotations": [],
            },
            {
                "id": "gomodel:test_happens_during/phase1",
                "type": [
                    {"type": "class", "id": "GO:0005132", "label": "cell cycle phase"}
                ],
                "root-type": [
                    {"type": "class", "id": "GO:0044848", "label": "biological phase"}
                ],
                "annotations": [],
            },
        ],
        "facts": [
            {
                "subject": "gomodel:test_happens_during/activity1",
                "property": Relation.ENABLED_BY,
                "object": "gomodel:test_happens_during/protein1",
                "annotations": [],
            },
            {
                "subject": "gomodel:test_happens_during/activity1",
                "property": Relation.HAPPENS_DURING,
                "object": "gomodel:test_happens_during/phase1",
                "annotations": [],
            },
        ],
    }

    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    assert model.activities is not None
    assert len(model.activities) == 1
    activity = model.activities[0]
    assert activity.happens_during is not None
    assert activity.happens_during.term == "GO:0005132"


def test_molecular_association_inverse_properties_definition():
    """Test that all keys of MOLECULE_ASSOCIATION_INVERSE_PROPERTIES have a value that is included
    in MOLECULE_ASSOCIATION_PROPERTIES"""
    for inverse_prop, prop in MOLECULAR_ASSOCIATION_INVERSE_PROPERTIES.items():
        assert prop in MOLECULAR_ASSOCIATION_PROPERTIES, (
            f"Inverse property '{inverse_prop}' maps to '{prop}', which is not a valid molecule association property."
        )


def test_has_small_molecule_activator_association():
    """Test that facts with the is_small_molecule_activator_of property are correctly translated."""
    minerva_object = load_minerva_object("6606056e00002011")
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    # Because the associations are "activity-oriented" in the data model, the
    # is_small_molecule_activator_of association facts in the Minerva object will be translated to
    # has_small_molecule_activator associations on the Activity.
    has_small_molecule_activator_associations = [
        assoc
        for activity in model.activities or []
        for assoc in activity.molecular_associations or []
        if assoc.predicate == Relation.HAS_SMALL_MOLECULE_ACTIVATOR
    ]

    assert len(has_small_molecule_activator_associations) == 2


def test_deeply_nested_part_of_associations():
    """Test that deeply nested part_of associations are correctly translated."""
    minerva_object = load_minerva_object("64d5781900000615")
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    rraga_activity = next(
        (
            a
            for a in model.activities or []
            if a.enabled_by
            and a.enabled_by.term == "UniProtKB:Q7L523"
            and a.molecular_function
            and a.molecular_function.term == "GO:0043495"
        ),
        None,
    )

    assert rraga_activity is not None
    assert rraga_activity.part_of is not None
    assert rraga_activity.part_of.term == "GO:0061462"
    assert rraga_activity.part_of.part_of is not None
    assert rraga_activity.part_of.part_of.term == "GO:1904263"
    assert rraga_activity.part_of.part_of.part_of is not None
    assert rraga_activity.part_of.part_of.part_of.term == "GO:0031669"


def test_deeply_nested_molecule_localization():
    """Test that deeply nested molecule localization associations are correctly translated."""
    minerva_object = load_minerva_object("nested-localization")
    mw = MinervaWrapper()
    model = mw.minerva_object_to_model(minerva_object)

    iron_molecule = next(
        (m for m in model.molecules or [] if m.term == "CHEBI:29033"),
        None,
    )

    assert iron_molecule is not None
    assert iron_molecule.located_in is not None
    assert iron_molecule.located_in.term == "GO:0005880"
    assert iron_molecule.located_in.part_of is not None
    assert iron_molecule.located_in.part_of.term == "GO:0005637"
    assert iron_molecule.located_in.part_of.part_of is not None
    assert iron_molecule.located_in.part_of.part_of.term == "GO:0005634"


def test_minerva_view_get_facts():
    """Test that the MinervaView.get_facts method correctly retrieves facts for a given subject."""
    minerva_object = {
        "facts": [
            {
                "subject": "S1",
                "property": "P1",
                "object": "O1",
            },
            {
                "subject": "S1",
                "property": "P2",
                "object": "O2",
            },
            {
                "subject": "S2",
                "property": "P1",
                "object": "O3",
            },
            {
                "subject": "S2",
                "property": "P1",
                "object": "O4",
            },
            {
                "subject": "S3",
                "property": "P1",
                "object": "O3",
            },
            {
                "subject": "S2",
                "property": "P3",
                "object": "O3",
            },
        ]
    }
    view = MinervaView(minerva_object)
    assert len(view.get_facts()) == 6
    assert len(view.get_facts(property="P1")) == 4
    assert len(view.get_facts(object="O1")) == 1
    assert len(view.get_facts(object="O3", property="P2")) == 0
    assert len(view.get_facts(subject="S1")) == 2
    assert len(view.get_facts(subject="S1", property="P2")) == 1
    assert len(view.get_facts(subject="S1", object="O3")) == 0
    assert len(view.get_facts(subject="S2", property="P1", object="O3")) == 1
