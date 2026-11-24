# Derived from:
# https://github.com/geneontology/web-components/blob/5d87e593121eafe6ac4690fa4591f88aa5a03fd8/packages/web-components/src/globals/%40noctua.form/data/taxon-dataset.json
import logging
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
    EnabledByGeneProductSetAssociation,
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

logger = logging.getLogger(__name__)

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

# If one activity has a molecular association with predicate P1 and molecule M,
# and another activity has a molecular association with predicate P2 and the same
# molecule M, then we can infer a causal relationship between the two activities. The
# following set of predicate pairs defines which combinations of predicates can be
# used to infer such causal relationships.
IMPLICIT_CAUSAL_ASSOCIATION_CHAINS: frozenset[tuple[Relation, Relation]] = frozenset(
    {
        (Relation.HAS_OUTPUT, Relation.HAS_INPUT),
        (Relation.HAS_OUTPUT, Relation.HAS_PRIMARY_INPUT),
        (Relation.HAS_OUTPUT, Relation.HAS_SMALL_MOLECULE_REGULATOR),
        (Relation.HAS_OUTPUT, Relation.HAS_SMALL_MOLECULE_ACTIVATOR),
        (Relation.HAS_OUTPUT, Relation.HAS_SMALL_MOLECULE_INHIBITOR),
        (Relation.HAS_PRIMARY_OUTPUT, Relation.HAS_INPUT),
        (Relation.HAS_PRIMARY_OUTPUT, Relation.HAS_PRIMARY_INPUT),
        (
            Relation.HAS_PRIMARY_OUTPUT,
            Relation.HAS_SMALL_MOLECULE_REGULATOR,
        ),
        (
            Relation.HAS_PRIMARY_OUTPUT,
            Relation.HAS_SMALL_MOLECULE_ACTIVATOR,
        ),
        (
            Relation.HAS_PRIMARY_OUTPUT,
            Relation.HAS_SMALL_MOLECULE_INHIBITOR,
        ),
    }
)


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

        case EnabledByGeneProductAssociation():
            yield obj
            if obj.part_of:
                for assoc in obj.part_of:
                    yield from all_associations(assoc)

        case ProteinComplexHasPartAssociation():
            yield obj
            yield from obj.members or []
            if obj.part_of:
                for assoc in obj.part_of:
                    yield from all_associations(assoc)

        case EnabledByGeneProductSetAssociation():
            yield obj
            yield from obj.members

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


def enabled_by_gene_products(enabled_by: EnabledByAssociation) -> list[str]:
    """
    Get the gene products that (may) carry out an activity, given its enabled_by association.

    For a single gene product this is the gene product itself; for a protein complex it is the
    gene products that are parts of the complex (including members of sets that are parts of the
    complex); and for a set of substitutable gene products it is the members of the set.

    Example:
        >>> from gocam.datamodel import GeneProductSetMemberAssociation
        >>> erk_set = EnabledByGeneProductSetAssociation(
        ...     term="CHEBI:33695",
        ...     members=[
        ...         GeneProductSetMemberAssociation(term="UniProtKB:P27361"),  # ERK1
        ...         GeneProductSetMemberAssociation(term="UniProtKB:P28482"),  # ERK2
        ...     ],
        ... )
        >>> enabled_by_gene_products(erk_set)
        ['UniProtKB:P27361', 'UniProtKB:P28482']
        >>> enabled_by_gene_products(EnabledByGeneProductAssociation(term="UniProtKB:P27361"))
        ['UniProtKB:P27361']
        >>> complex_with_set = EnabledByProteinComplexAssociation(
        ...     term="GO:0032991",
        ...     has_part=[
        ...         ProteinComplexHasPartAssociation(term="UniProtKB:Q02750"),  # MEK1
        ...         ProteinComplexHasPartAssociation(term="CHEBI:33695", members=erk_set.members),
        ...     ],
        ... )
        >>> enabled_by_gene_products(complex_with_set)
        ['UniProtKB:Q02750', 'UniProtKB:P27361', 'UniProtKB:P28482']

    Args:
        enabled_by: The enabled_by association of an activity.

    Returns:
        A list of gene product IDs.
    """
    match enabled_by:
        case EnabledByGeneProductSetAssociation():
            return [m.term for m in enabled_by.members if m.term]
        case EnabledByProteinComplexAssociation():
            gene_products = []
            for part in enabled_by.has_part or []:
                if part.members:
                    gene_products.extend(m.term for m in part.members if m.term)
                elif part.term:
                    gene_products.append(part.term)
            return gene_products
        case _:
            return [enabled_by.term] if enabled_by.term else []


@overload
def all_provenance(obj: Model) -> Iterator[ProvenanceInfo]: ...
@overload
def all_provenance(obj: Association) -> Iterator[ProvenanceInfo]: ...
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

    Nodes represent enabled activities. Edges represent explicit causal
    associations or molecule-mediated predicate chains declared in
    IMPLICIT_CAUSAL_ASSOCIATION_CHAINS.

    Each edge records how it arose, as two independent boolean attributes:

        ``explicit``  — asserted by a causal association on the activity.
        ``inferred``  — implied by two activities sharing a molecule across a
                        predicate chain.

    They are independent because a pair of activities can be connected both
    ways at once, so callers wanting only the implied connections must filter
    on ``inferred`` rather than subtracting the explicit ones.

    Args:
        model: The model to convert.

    Returns:
        A directed graph representing causal connectivity in the model.
    """
    graph = nx.DiGraph()

    activities_by_molecule_and_predicate: defaultdict[
        tuple[str, Relation], set[str]
    ] = defaultdict(set)
    for activity in model.activities or []:
        for association in activity.molecular_associations or []:
            if association.molecule is None:
                continue
            try:
                predicate = Relation(association.predicate)
            except ValueError:
                logger.warning("Unknown predicate: %s", association.predicate)
                continue
            activities_by_molecule_and_predicate[(association.molecule, predicate)].add(
                activity.id
            )

    enabled_activity_ids = {
        a.id for a in model.activities or [] if a.enabled_by is not None
    }

    for activity in model.activities or []:
        if activity.id not in enabled_activity_ids:
            continue

        graph.add_node(activity.id)

        for causal_assoc in activity.causal_associations or []:
            downstream_activity_id = causal_assoc.downstream_activity
            if downstream_activity_id in enabled_activity_ids:
                graph.add_edge(activity.id, downstream_activity_id, explicit=True)

        for association in activity.molecular_associations or []:
            if association.molecule is None:
                continue
            for (
                upstream_predicate,
                downstream_predicate,
            ) in IMPLICIT_CAUSAL_ASSOCIATION_CHAINS:
                if association.predicate != upstream_predicate:
                    continue
                downstream_activity_ids = activities_by_molecule_and_predicate.get(
                    (association.molecule, downstream_predicate), set()
                )
                for downstream_activity_id in downstream_activity_ids:
                    if (
                        downstream_activity_id != activity.id
                        and downstream_activity_id in enabled_activity_ids
                    ):
                        graph.add_edge(
                            activity.id, downstream_activity_id, inferred=True
                        )

    return graph
