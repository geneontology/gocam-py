import pytest

from gocam.datamodel import (
    Activity,
    CausalAssociation,
    EnabledByGeneProductAssociation,
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
    assert graph.edges["upstream", "downstream"] == {}


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
