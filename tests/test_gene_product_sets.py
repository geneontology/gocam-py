"""Tests for activities enabled by sets of substitutable gene products (e.g. paralogy groups).

See https://github.com/geneontology/gocam-py/issues/68 and
https://github.com/geneontology/pathways2GO/issues/352.

The ``minerva-R-HSA-*.json`` fixtures are Reactome test models produced by pathways2GO with set
support (geneontology/noctua-models PR #287, ``dev-2025-06-04`` branch). Those models are not
currently served by any Minerva instance, so the fixtures were converted from the TTL files in
noctua-models to Minerva JSON; class labels and root types (which Minerva computes from the
ontology) were filled in from the production Minerva JSON for the same models.
"""

import doctest
import json

import pytest
import yaml
from pydantic import ValidationError

import gocam.utils
from gocam.datamodel import (
    EnabledByGeneProductAssociation,
    EnabledByGeneProductSetAssociation,
    EnabledByProteinComplexAssociation,
    GeneProductSetMemberAssociation,
    Model,
)
from gocam.translation.minerva_wrapper import MinervaWrapper
from gocam.translation.result import WarningType
from gocam.utils import all_associations, enabled_by_gene_products
from gocam.vocabulary import INFORMATION_BIOMACROMOLECULE, Relation
from tests import INPUT_DIR

ERK1 = "UniProtKB:P27361"
ERK2 = "UniProtKB:P28482"
MEK1 = "UniProtKB:Q02750"


def load_minerva_object(id: str) -> dict:
    with open(INPUT_DIR / f"minerva-{id}.json") as f:
        return json.load(f)


@pytest.fixture
def erk_model() -> Model:
    with open(INPUT_DIR / "Model-erk-paralog-set.yaml") as f:
        return Model.model_validate(yaml.safe_load(f))


@pytest.fixture
def reactome_set_model() -> Model:
    """Reactome model in which three of the four activities are enabled by sets."""
    return MinervaWrapper.minerva_object_to_model(load_minerva_object("R-HSA-9634600"))


def _activity(model: Model, activity_id: str):
    return next(a for a in model.activities or [] if a.id == activity_id)


# -- data model ---------------------------------------------------------------------------------


def test_set_from_yaml(erk_model):
    """The type designator selects the set association when loading from YAML/JSON."""
    eb = _activity(erk_model, "gomodel:example_paralogy_model/erk_activity").enabled_by
    assert isinstance(eb, EnabledByGeneProductSetAssociation)
    assert eb.term == INFORMATION_BIOMACROMOLECULE
    assert eb.label == "ERK1/2"
    assert [m.term for m in eb.members] == [ERK1, ERK2]
    # set-level evidence applies to all members
    assert [e.reference for e in eb.evidence or []] == ["PMID:12345678"]


def test_set_json_round_trip(erk_model):
    round_tripped = Model.model_validate_json(erk_model.model_dump_json())
    assert round_tripped == erk_model


def test_set_requires_members():
    with pytest.raises(ValidationError):
        EnabledByGeneProductSetAssociation(term=INFORMATION_BIOMACROMOLECULE)  # ty:ignore[missing-argument]


def test_all_associations_includes_set_members(erk_model):
    associations = list(all_associations(erk_model))
    member_terms = [
        a.term for a in associations if isinstance(a, GeneProductSetMemberAssociation)
    ]
    assert member_terms == [ERK1, ERK2]


def test_utils_doctests():
    results = doctest.testmod(gocam.utils)
    assert results.attempted > 0
    assert results.failed == 0


# -- Minerva translation ------------------------------------------------------------------------


def _minerva_set_model(set_facts_evidence: bool = True) -> dict:
    """A minimal hand-built Minerva model: an MF enabled by a set of ERK1 and ERK2."""

    def ev(local_id, eco, source):
        return {
            "id": f"gomodel:erk/{local_id}",
            "type": [{"type": "class", "id": eco, "label": eco}],
            "root-type": [{"type": "class", "id": "ECO:0000000", "label": "evidence"}],
            "annotations": [{"key": "source", "value": source}],
        }

    gp_root = [
        {"type": "class", "id": "CHEBI:33695", "label": "information biomacromolecule"},
        {"type": "class", "id": "CHEBI:24431", "label": "chemical entity"},
    ]
    individuals = [
        {
            "id": "gomodel:erk/mf",
            "type": [
                {"type": "class", "id": "GO:0004707", "label": "MAP kinase activity"}
            ],
            "root-type": [
                {"type": "class", "id": "GO:0003674", "label": "molecular_function"}
            ],
        },
        {
            "id": "gomodel:erk/set",
            "type": [
                {
                    "type": "class",
                    "id": "CHEBI:33695",
                    "label": "information biomacromolecule",
                }
            ],
            "root-type": gp_root,
            "annotations": [{"key": "rdfs:label", "value": "ERK1/2"}],
        },
        {
            "id": "gomodel:erk/erk1",
            "type": [{"type": "class", "id": ERK1, "label": "MAPK3 Hsap"}],
            "root-type": gp_root,
        },
        {
            "id": "gomodel:erk/erk2",
            "type": [{"type": "class", "id": ERK2, "label": "MAPK1 Hsap"}],
            "root-type": gp_root,
        },
        ev("ev_set", "ECO:0000304", "PMID:12345678"),
        ev("ev_erk1", "ECO:0000250", "PMID:1"),
    ]

    def evidence(local_id):
        return [
            {"key": "evidence", "value": f"gomodel:erk/{local_id}", "value-type": "IRI"}
        ]

    facts = [
        {
            "subject": "gomodel:erk/mf",
            "property": Relation.ENABLED_BY.value,
            "object": "gomodel:erk/set",
            "annotations": evidence("ev_set"),
        },
        {
            "subject": "gomodel:erk/set",
            "property": Relation.HAS_SUBSTITUTABLE_ENTITY.value,
            "object": "gomodel:erk/erk1",
            "annotations": evidence("ev_erk1"),
        },
        {
            "subject": "gomodel:erk/set",
            "property": Relation.HAS_SUBSTITUTABLE_ENTITY.value,
            "object": "gomodel:erk/erk2",
            "annotations": [],
        },
    ]
    return {
        "id": "gomodel:erk",
        "individuals": individuals,
        "facts": facts,
        "annotations": [{"key": "title", "value": "ERK set"}],
    }


def test_minerva_set_minimal():
    result = MinervaWrapper.translate(_minerva_set_model())
    assert result.warnings == []
    [activity] = result.result.activities or []
    eb = activity.enabled_by
    assert isinstance(eb, EnabledByGeneProductSetAssociation)
    assert eb.term == INFORMATION_BIOMACROMOLECULE
    assert eb.label == "ERK1/2"
    # evidence on the enabled_by edge is the evidence for the set as a whole
    assert [(e.term, e.reference) for e in eb.evidence or []] == [
        ("ECO:0000304", "PMID:12345678")
    ]
    # evidence on a has substitutable entity edge stays with that membership assertion
    erk1, erk2 = eb.members
    assert erk1.term == ERK1
    assert [e.reference for e in erk1.evidence or []] == ["PMID:1"]
    assert erk2.term == ERK2
    assert erk2.evidence == []


def test_minerva_reactome_sets(reactome_set_model):
    activities = reactome_set_model.activities or []
    assert len(activities) == 4
    set_activities = [
        a
        for a in activities
        if isinstance(a.enabled_by, EnabledByGeneProductSetAssociation)
    ]
    assert len(set_activities) == 3

    eb = _activity(reactome_set_model, "gomodel:R-HSA-163773").enabled_by
    assert isinstance(eb, EnabledByGeneProductSetAssociation)
    assert eb.label == "PRKACA,PRKACB,PRKACG"
    assert sorted(m.term for m in eb.members) == [
        "UniProtKB:P17612",
        "UniProtKB:P22612",
        "UniProtKB:P22694",
    ]
    assert {e.reference for e in eb.evidence or []} == {
        "Reactome:R-HSA-163773",
        "Reactome:R-HSA-9634600",
    }

    # an ordinary gene product enabler in the same model is unaffected
    eb = _activity(reactome_set_model, "gomodel:R-HSA-70262").enabled_by
    assert isinstance(eb, EnabledByGeneProductAssociation)
    assert eb.term == "UniProtKB:P16118"


def test_minerva_reactome_sets_no_warnings():
    """All has substitutable entity facts are consumed by the translation."""
    result = MinervaWrapper.translate(load_minerva_object("R-HSA-9634600"))
    assert result.warnings == []


def test_minerva_reactome_complex_with_set_parts():
    """Sets can also be parts of a complex that enables an activity."""
    result = MinervaWrapper.translate(load_minerva_object("R-HSA-1855229"))
    eb = _activity(result.result, "gomodel:R-HSA-1855188").enabled_by
    assert isinstance(eb, EnabledByProteinComplexAssociation)
    set_parts = [p for p in eb.has_part or [] if p.members]
    assert len(set_parts) == 2
    assert all(p.term == INFORMATION_BIOMACROMOLECULE for p in set_parts)
    assert sorted(m.term for p in set_parts for m in p.members or []) == [
        "UniProtKB:A8CG34",
        "UniProtKB:Isoform_P52948-3",
        "UniProtKB:Isoform_P52948-4",
        "UniProtKB:Q96HA1",
    ]
    assert "UniProtKB:A8CG34" in enabled_by_gene_products(eb)
    assert INFORMATION_BIOMACROMOLECULE not in enabled_by_gene_products(eb)
    # The remaining unhandled facts are about complexes nested inside complexes, which are
    # not yet representable, not about sets directly used as parts of the enabler.
    for warning in result.warnings:
        assert warning.type == WarningType.UNHANDLED_FACT


# -- other translations -------------------------------------------------------------------------


def test_tbox_set(erk_model):
    pytest.importorskip("pyhornedowl")
    from gocam.translation.tbox_translator import TBoxTranslator

    translator = TBoxTranslator()
    translator.load_models([erk_model])
    axioms = [str(a) for a in translator.ontology.get_axioms()]
    substitutable = [a for a in axioms if "RO_0019003" in a]
    assert len(substitutable) == 2
    assert any("P27361" in a for a in substitutable)
    assert any("P28482" in a for a in substitutable)


def test_cx2_set(erk_model):
    pytest.importorskip("ndex2")
    from gocam.translation.cx2 import model_to_cx2

    cx2 = model_to_cx2(erk_model)
    nodes = next(aspect["nodes"] for aspect in cx2 if "nodes" in aspect)
    set_node = next(n for n in nodes if n["v"]["type"] == "gene set")
    assert set_node["v"]["name"] == "ERK1/2"
    assert set_node["v"]["member"] == ["MAPK3", "MAPK1"]
    assert set_node["v"]["represents"] == INFORMATION_BIOMACROMOLECULE


def test_networkx_set(erk_model):
    from gocam.translation.networkx.model_network_translator import (
        ModelNetworkTranslator,
    )

    graph = ModelNetworkTranslator().translate_models([erk_model])
    # each member of the set becomes a gene node, rather than one node for CHEBI:33695
    assert set(graph.nodes) == {MEK1, ERK1, ERK2}
    assert set(graph.edges) == {(MEK1, ERK1), (MEK1, ERK2)}
    assert graph.nodes[ERK1]["label"] == "MAPK3 Hsap"


@pytest.mark.parametrize(
    "enabled_by,expected",
    [
        (EnabledByGeneProductAssociation(term=MEK1), [MEK1]),
        (
            EnabledByGeneProductSetAssociation(
                term=INFORMATION_BIOMACROMOLECULE,
                members=[
                    GeneProductSetMemberAssociation(term=ERK1),
                    GeneProductSetMemberAssociation(term=ERK2),
                ],
            ),
            [ERK1, ERK2],
        ),
    ],
)
def test_enabled_by_gene_products(enabled_by, expected):
    assert enabled_by_gene_products(enabled_by) == expected
