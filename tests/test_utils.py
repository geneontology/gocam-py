import pytest

from gocam.datamodel import (
    Activity,
    CausalAssociation,
    EnabledByGeneProductAssociation,
    EvidenceItem,
    Model,
    MoleculeAssociation,
)
from gocam.utils import model_to_digraph
from gocam.vocabulary import Relation

EXPECTED_IMPLICIT_CAUSAL_ASSOCIATION_CHAINS = (
    (Relation.HAS_OUTPUT, Relation.HAS_INPUT),
    (Relation.HAS_OUTPUT, Relation.HAS_PRIMARY_INPUT),
    (Relation.HAS_OUTPUT, Relation.HAS_SMALL_MOLECULE_REGULATOR),
    (Relation.HAS_OUTPUT, Relation.HAS_SMALL_MOLECULE_ACTIVATOR),
    (Relation.HAS_OUTPUT, Relation.HAS_SMALL_MOLECULE_INHIBITOR),
    (Relation.HAS_PRIMARY_OUTPUT, Relation.HAS_INPUT),
    (Relation.HAS_PRIMARY_OUTPUT, Relation.HAS_PRIMARY_INPUT),
    (Relation.HAS_PRIMARY_OUTPUT, Relation.HAS_SMALL_MOLECULE_REGULATOR),
    (Relation.HAS_PRIMARY_OUTPUT, Relation.HAS_SMALL_MOLECULE_ACTIVATOR),
    (Relation.HAS_PRIMARY_OUTPUT, Relation.HAS_SMALL_MOLECULE_INHIBITOR),
)


def activity_with_molecular_associations(
    activity_id: str,
    *associations: tuple[Relation, str],
) -> Activity:
    return Activity(
        id=activity_id,
        enabled_by=EnabledByGeneProductAssociation(term=f"UniProtKB:{activity_id}"),
        molecular_associations=[
            MoleculeAssociation(predicate=predicate, molecule=molecule)
            for predicate, molecule in associations
        ],
    )


@pytest.mark.parametrize(
    ("upstream_predicate", "downstream_predicate"),
    EXPECTED_IMPLICIT_CAUSAL_ASSOCIATION_CHAINS,
)
def test_declared_implicit_causal_chain_adds_edge(
    upstream_predicate: Relation,
    downstream_predicate: Relation,
):
    model = Model(
        id="gomodel:test",
        title="Implicit causal chain test",
        activities=[
            activity_with_molecular_associations(
                "upstream", (upstream_predicate, "molecule")
            ),
            activity_with_molecular_associations(
                "downstream", (downstream_predicate, "molecule")
            ),
        ],
    )

    graph = model_to_digraph(model)

    assert set(graph.nodes) == {"upstream", "downstream"}
    assert set(graph.edges) == {("upstream", "downstream")}
    assert graph.edges["upstream", "downstream"] == {"inferred": True}


@pytest.mark.parametrize(
    ("upstream_predicate", "downstream_predicate"),
    [
        (Relation.HAS_INPUT, Relation.HAS_INPUT),
        (Relation.HAS_OUTPUT, Relation.HAS_OUTPUT),
        (
            Relation.HAS_SMALL_MOLECULE_ACTIVATOR,
            Relation.HAS_INPUT,
        ),
    ],
)
def test_undeclared_molecular_association_pair_does_not_add_edge(
    upstream_predicate: Relation,
    downstream_predicate: Relation,
):
    model = Model(
        id="gomodel:test",
        title="Undeclared chain test",
        activities=[
            activity_with_molecular_associations(
                "first", (upstream_predicate, "molecule")
            ),
            activity_with_molecular_associations(
                "second", (downstream_predicate, "molecule")
            ),
        ],
    )

    assert model_to_digraph(model).number_of_edges() == 0


def test_molecular_associations_to_different_molecules_do_not_add_edge():
    model = Model(
        id="gomodel:test",
        title="Different molecule test",
        activities=[
            activity_with_molecular_associations(
                "upstream", (Relation.HAS_OUTPUT, "molecule-1")
            ),
            activity_with_molecular_associations(
                "downstream", (Relation.HAS_INPUT, "molecule-2")
            ),
        ],
    )

    assert model_to_digraph(model).number_of_edges() == 0


def test_molecular_association_without_molecule_does_not_add_edge():
    upstream = activity_with_molecular_associations("upstream")
    upstream.molecular_associations = [
        MoleculeAssociation(predicate=Relation.HAS_OUTPUT)
    ]
    model = Model(
        id="gomodel:test",
        title="Missing molecule test",
        activities=[
            upstream,
            activity_with_molecular_associations(
                "downstream", (Relation.HAS_INPUT, "molecule")
            ),
        ],
    )

    assert model_to_digraph(model).number_of_edges() == 0


def test_implicit_causal_chain_does_not_add_self_edge():
    model = Model(
        id="gomodel:test",
        title="Implicit self-edge test",
        activities=[
            activity_with_molecular_associations(
                "activity",
                (Relation.HAS_OUTPUT, "molecule"),
                (Relation.HAS_INPUT, "molecule"),
            )
        ],
    )

    graph = model_to_digraph(model)

    assert set(graph.nodes) == {"activity"}
    assert graph.number_of_edges() == 0


def test_duplicate_implicit_chains_collapse_to_one_edge():
    model = Model(
        id="gomodel:test",
        title="Duplicate implicit edge test",
        activities=[
            activity_with_molecular_associations(
                "upstream",
                (Relation.HAS_OUTPUT, "molecule"),
                (Relation.HAS_OUTPUT, "molecule"),
            ),
            activity_with_molecular_associations(
                "downstream",
                (Relation.HAS_INPUT, "molecule"),
                (Relation.HAS_INPUT, "molecule"),
            ),
        ],
    )

    graph = model_to_digraph(model)

    assert set(graph.edges) == {("upstream", "downstream")}
    assert graph.number_of_edges() == 1


def test_explicit_causal_association_still_adds_edge():
    upstream = activity_with_molecular_associations("upstream")
    upstream.causal_associations = [
        CausalAssociation(
            predicate=Relation.CAUSALLY_UPSTREAM_OF,
            downstream_activity="downstream",
        )
    ]
    model = Model(
        id="gomodel:test",
        title="Explicit causal edge test",
        activities=[
            upstream,
            activity_with_molecular_associations("downstream"),
        ],
    )

    graph = model_to_digraph(model)

    assert set(graph.edges) == {("upstream", "downstream")}
    assert graph.edges["upstream", "downstream"] == {"explicit": True}


def test_optional_support_retains_predicates_intermediates_and_evidence():
    upstream = activity_with_molecular_associations(
        "upstream", (Relation.HAS_OUTPUT, "m1"), (Relation.HAS_PRIMARY_OUTPUT, "m2")
    )
    downstream = activity_with_molecular_associations(
        "downstream",
        (Relation.HAS_INPUT, "m1"),
        (Relation.HAS_SMALL_MOLECULE_INHIBITOR, "m2"),
    )
    upstream.causal_associations = [
        CausalAssociation(
            predicate=Relation.CAUSALLY_UPSTREAM_OF,
            downstream_activity="downstream",
            evidence=[EvidenceItem(reference="PMID:1")],
        ),
        CausalAssociation(
            predicate="RO:0002413",
            downstream_activity="downstream",
            evidence=[EvidenceItem(reference="PMID:2")],
        ),
    ]
    assert upstream.molecular_associations is not None
    assert downstream.molecular_associations is not None
    upstream.molecular_associations[0].evidence = [EvidenceItem(reference="PMID:3")]
    downstream.molecular_associations[0].evidence = [EvidenceItem(reference="PMID:4")]
    model = Model(
        id="gomodel:test",
        title="Supporting associations",
        activities=[upstream, downstream],
    )
    before = model.model_dump()
    graph = model_to_digraph(model, include_support=True)
    edge = graph.edges["upstream", "downstream"]
    assert edge["explicit"] and edge["inferred"]
    assert edge["causal_associations"] == upstream.causal_associations
    chains = {chain["molecule"]: chain for chain in edge["molecule_chains"]}
    assert set(chains) == {"m1", "m2"}
    assert chains["m1"]["upstream"] is upstream.molecular_associations[0]
    assert chains["m1"]["downstream"] is downstream.molecular_associations[0]
    assert chains["m1"]["upstream"].evidence[0].reference == "PMID:3"
    assert chains["m1"]["downstream"].evidence[0].reference == "PMID:4"
    assert chains["m2"]["downstream"].predicate == Relation.HAS_SMALL_MOLECULE_INHIBITOR
    assert model.model_dump() == before
    assert model_to_digraph(model).edges["upstream", "downstream"] == {
        "explicit": True,
        "inferred": True,
    }


def test_support_groups_identical_chains_but_keeps_different_evidence():
    upstream = activity_with_molecular_associations(
        "upstream", (Relation.HAS_OUTPUT, "m")
    )
    assert upstream.molecular_associations is not None
    output = upstream.molecular_associations[0]
    upstream.molecular_associations.extend(
        [output.model_copy(deep=True), output.model_copy(deep=True)]
    )
    upstream.molecular_associations[-1].evidence = [EvidenceItem(reference="PMID:5")]
    downstream = activity_with_molecular_associations(
        "downstream", (Relation.HAS_INPUT, "m")
    )
    model = Model(
        id="gomodel:test", title="Duplicate support", activities=[upstream, downstream]
    )
    graph = model_to_digraph(model, include_support=True)
    assert len(graph.edges["upstream", "downstream"]["molecule_chains"]) == 2
