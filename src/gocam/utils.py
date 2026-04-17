# Derived from:
# https://github.com/geneontology/web-components/blob/5d87e593121eafe6ac4690fa4591f88aa5a03fd8/packages/web-components/src/globals/%40noctua.form/data/taxon-dataset.json
from collections import defaultdict
from typing import Any, Iterator, overload

import networkx as nx

from gocam.datamodel import (
    Activity,
    Association,
    BiologicalProcessAssociation,
    CellTypeAssociation,
    CellularAnatomicalEntityAssociation,
    EnabledByAssociation,
    EnabledByGeneProductAssociation,
    EnabledByProteinComplexAssociation,
    EvidenceItem,
    GrossAnatomyAssociation,
    Model,
    MoleculeAssociation,
    MoleculeNode,
    PartOfProteinComplexAssociation,
    ProteinComplexHasPartAssociation,
    ProvenanceInfo,
)
from gocam.vocabulary import Relation

SPECIES_CODES = [
    "Atal",
    "Btau",
    "Cele",
    "Cfam",
    "Ddis",
    "Dmel",
    "Drer",
    "Ggal",
    "Hsap",
    "Mmus",
    "Pseudomonas",
    "Rnor",
    "Scer",
    "Sjap",
    "Solanaceae",
    "Spom",
    "Sscr",
    "Xenopus",
]


def remove_species_code_suffix(label: str) -> str:
    """Remove known species codes from the end of a label.

    If a label ends with one of the known species codes, remove it and any trailing whitespace.
    Otherwise, return the label unchanged.

    :param label: The label to process.
    :return: The processed label.
    """
    for code in SPECIES_CODES:
        label = label.removesuffix(code).strip()
    return label


def all_activity_inputs(activity: Activity) -> list[MoleculeAssociation]:
    """Get all inputs for an activity, including primary input if present.

    Args:
        activity: The activity to get inputs for.

    Returns:
        List of all molecule associations that are inputs to the activity.
    """
    return [
        ma
        for ma in activity.molecular_associations or []
        if ma.predicate == Relation.HAS_INPUT
        or ma.predicate == Relation.HAS_PRIMARY_INPUT
    ]


def all_activity_outputs(activity: Activity) -> list[MoleculeAssociation]:
    """Get all outputs for an activity, including primary output if present.

    Args:
        activity: The activity to get outputs for.

    Returns:
        List of all molecule associations that are outputs of the activity.
    """
    return [
        ma
        for ma in activity.molecular_associations or []
        if ma.predicate == Relation.HAS_OUTPUT
        or ma.predicate == Relation.HAS_PRIMARY_OUTPUT
    ]


@overload
def all_associations(obj: Model) -> Iterator[Association]: ...
@overload
def all_associations(obj: Activity) -> Iterator[Association]: ...
@overload
def all_associations(obj: MoleculeNode) -> Iterator[Association]: ...
@overload
def all_associations(obj: EnabledByAssociation) -> Iterator[Association]: ...
@overload
def all_associations(obj: BiologicalProcessAssociation) -> Iterator[Association]: ...
@overload
def all_associations(
    obj: CellularAnatomicalEntityAssociation,
) -> Iterator[Association]: ...
@overload
def all_associations(obj: CellTypeAssociation) -> Iterator[Association]: ...
@overload
def all_associations(obj: GrossAnatomyAssociation) -> Iterator[Association]: ...
@overload
def all_associations(obj: Association) -> Iterator[Association]: ...
def all_associations(obj: Any) -> Iterator[Association]:
    """
    Extract all Association objects from a given object.

    Args:
        obj: The object to extract associations from. Can be a Model, Activity, MoleculeNode, or
             any Association type.

    Yields:
        An Association object extracted from the input object or its nested associations.
    """
    match obj:
        case Model():
            for activity in obj.activities or []:
                yield from all_associations(activity)
            for molecule in obj.molecules or []:
                yield from all_associations(molecule)

        case Activity():
            if obj.enabled_by:
                yield from all_associations(obj.enabled_by)
            if obj.molecular_function:
                yield from all_associations(obj.molecular_function)
            if obj.part_of:
                yield from all_associations(obj.part_of)
            if obj.occurs_in:
                yield from all_associations(obj.occurs_in)
            if obj.happens_during:
                yield from all_associations(obj.happens_during)
            for molecule_association in obj.molecular_associations or []:
                yield from all_associations(molecule_association)
            for causal_association in obj.causal_associations or []:
                yield from all_associations(causal_association)

        case MoleculeNode():
            if obj.located_in:
                yield from all_associations(obj.located_in)

        case EnabledByGeneProductAssociation() | ProteinComplexHasPartAssociation():
            yield obj
            if obj.part_of:
                for assoc in obj.part_of:
                    yield from all_associations(assoc)

        case EnabledByProteinComplexAssociation() | PartOfProteinComplexAssociation():
            yield obj
            if obj.has_part:
                for assoc in obj.has_part:
                    yield from all_associations(assoc)

        case BiologicalProcessAssociation():
            yield obj
            if obj.happens_during:
                yield from all_associations(obj.happens_during)
            if obj.part_of:
                yield from all_associations(obj.part_of)

        case (
            CellularAnatomicalEntityAssociation()
            | CellTypeAssociation()
            | GrossAnatomyAssociation()
        ):
            yield obj
            if obj.part_of:
                yield from all_associations(obj.part_of)

        case Association():
            yield obj

        case _:
            raise ValueError(f"Unsupported object type: {type(obj)}")


@overload
def all_provenance(model: Model) -> Iterator[ProvenanceInfo]: ...
@overload
def all_provenance(association: Association) -> Iterator[ProvenanceInfo]: ...
def all_provenance(obj) -> Iterator[ProvenanceInfo]:
    """
    Extract all ProvenanceInfo objects from a given object

    Args:
        obj: The object to extract provenance information from. Can be a Model or any Association
             type.

    Yields:
        A ProvenanceInfo object extracted from the input object or its nested associations.
    """
    match obj:
        case Model():
            if obj.provenances:
                yield from obj.provenances
            for activity in obj.activities or []:
                if activity.provenances:
                    yield from activity.provenances
            for association in all_associations(obj):
                yield from all_provenance(association)

        case Association():
            if obj.provenances:
                yield from obj.provenances
            if obj.evidence:
                for evidence in obj.evidence:
                    if evidence.provenances:
                        yield from evidence.provenances

        case _:
            raise ValueError(f"Unsupported object type: {type(obj)}")


def all_evidence(model: Model) -> Iterator[EvidenceItem]:
    """
    Extract all EvidenceItem objects from a Model.

    Args:
        model: The GO-CAM model to extract evidence from.

    Yields:
        An EvidenceItem object extracted from the model's associations.
    """
    for association in all_associations(model):
        if association.evidence:
            yield from association.evidence


def model_to_digraph(model: Model) -> nx.DiGraph:
    """Convert a Model to a NetworkX directed graph.

    Nodes are added for activities with an 'enabled_by' property. Edges are added based on:
    - Causal associations between activities.
    - Outputs of one activity serving as inputs to another activity.

    Args:
        model: The model to convert.

    Returns:
        A directed graph representation of the model.
    """
    graph = nx.DiGraph()

    activities_by_input: dict[str, list[str]] = defaultdict(list)
    for activity in model.activities or []:
        for input_ in all_activity_inputs(activity):
            if input_.molecule:
                activities_by_input[input_.molecule].append(activity.id)

    for activity in model.activities or []:
        if activity.enabled_by is None:
            continue

        graph.add_node(activity.id)

        for causal_assoc in activity.causal_associations or []:
            downstream_activity_id = causal_assoc.downstream_activity
            if downstream_activity_id:
                graph.add_edge(activity.id, downstream_activity_id)

        for output in all_activity_outputs(activity):
            if output.molecule:
                for downstream_activity_id in activities_by_input.get(
                    output.molecule, []
                ):
                    if downstream_activity_id != activity.id:
                        graph.add_edge(activity.id, downstream_activity_id)

    return graph
