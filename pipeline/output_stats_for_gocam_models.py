"""
Output statistics information about GO-CAM models.

This module processes GO-CAM model JSON files and calculates a variety of
statistics at the per-model, per-contributor, per-provider, and aggregate
levels. Statistics include:

- Activity unit counts (total, enabled by gene product, enabled by protein complex)
- Gene product and protein complex enabler counts (total and unique)
- Causal relation counts
- Input/output molecule counts (chemical vs other, total vs unique)
- GO term counts (total and unique)
- Reference and PMID counts
- Inferred relation counts: An inferred relation is identified when an
  activity A produces chemical outputs (CHEBI terms) that are consumed as
  inputs by another activity B, and activity B also produces chemical outputs.
  This captures implicit metabolic pathway connections where one activity's
  product feeds into another activity's transformation. Each inferred relation
  records the pair of activity IDs [activity_A, activity_B].

Results are written as JSON files organized into subdirectories by model,
contributor (curator), and provider (group), along with aggregate summaries.
"""

import logging
import os
import re
from enum import Enum
from pathlib import Path
from typing import (
    Annotated,
    Collection,
    Dict,
    List,
    Literal,
    Set,
    TypeAlias,
)
from urllib.parse import urlparse

import typer
from rich import print
from rich.logging import RichHandler
from rich.progress import track
from rich.tree import Tree

from gocam.datamodel import (
    Association,
    ConfiguredBaseModel,
    EnabledByGeneProductAssociation,
    EnabledByProteinComplexAssociation,
    Model,
    QueryIndex,
)
from gocam.indexing.indexer import Indexer
from gocam.vocabulary.relation import Relation

app = typer.Typer()

logger = logging.getLogger(__name__)

INPUT_PREDICATES = {Relation.HAS_INPUT.value, Relation.HAS_PRIMARY_INPUT.value}
OUTPUT_PREDICATES = {Relation.HAS_OUTPUT.value, Relation.HAS_PRIMARY_OUTPUT.value}


class ResultType(str, Enum):
    SUCCESS = "Success"
    FILTERED = "Filtered"
    ERROR = "Error"


class FilterReason(str, Enum):
    """Reason a GO-CAM model was filtered out during processing."""

    NO_MODEL = "No model"


class ErrorReason(str, Enum):
    """Reason a GO-CAM model failed during processing."""

    READ_ERROR = "Read error"
    STATS_ERROR = "Error calculating statistics"
    WRITE_ERROR = "Write error"


class GocamStats(ConfiguredBaseModel):
    """Per-model or per-entity (contributor/provider) statistics for GO-CAM data.

    Numeric fields hold computed counts. Set and list fields are working
    collections used to accumulate data during processing; they feed the
    final numeric counts and are serialized for detailed output.
    """

    uri: str | None = None
    number_of_models: int = 0
    number_of_activity_units: int = 0
    number_of_activity_units_enabled_by_gene_product: int = 0
    number_of_activity_units_enabled_by_protein_complex: int = 0
    number_of_unique_gene_product_enablers: int = 0
    number_of_unique_protein_complex_genes: int = 0
    number_of_unique_gene_product_and_protein_complex_gene_enablers: int = 0
    number_of_genes: int = 0
    number_of_unique_references: int = 0
    number_of_unique_pmid: int = 0
    number_of_causal_relations: int = 0
    number_of_unique_causal_relations: int = 0
    number_of_inputs: int = 0
    number_of_chemical_inputs: int = 0
    number_of_other_inputs: int = 0
    number_of_unique_chemical_inputs: int = 0
    number_of_unique_other_inputs: int = 0
    number_of_outputs: int = 0
    number_of_chemical_outputs: int = 0
    number_of_other_outputs: int = 0
    number_of_unique_chemical_outputs: int = 0
    number_of_unique_other_outputs: int = 0
    number_of_go_terms: int = 0
    number_of_unique_go_terms: int = 0
    number_of_unique_inferred_relations: int = 0
    list_inferred_relations: List[List[str]] | None = []

    set_models: Set[str] = set()
    set_activities: Set[str] = set()
    set_enabled_by_gene_product: Set[str] = set()
    set_protein_complex_genes: Set[str] = set()
    list_enabled_by_gene_product: List[str] | None = []
    set_references: Set[str] = set()
    set_causal_relations: Set[str] = set()
    list_causal_relations: List[str] | None = []
    set_activity_unit_gene_product_enablers: Set[str] = set()
    set_activity_unit_protein_complex_enablers: Set[str] = set()
    set_protein_complex_in_activity_term: Set[str] = set()
    list_has_input_term: (
        List[str] | None
    ) = []  # Covers both has_input and has primary_input
    list_has_output_term: (
        List[str] | None
    ) = []  # Covers both has_output and has primary_output
    list_go_terms: List[str] | None = []


class ModelDetails(ConfiguredBaseModel):
    """Identifying metadata for a single processed GO-CAM model."""

    file_name: str = ""
    model_id: str = ""
    model_name: str = ""
    set_activities: Set[str] = set()


class AggregateInfo(ConfiguredBaseModel):
    """Aggregate statistics across all processed GO-CAM models, curators, or groups.

    Accumulates counts and working collections from individual GocamStats
    instances to produce cross-model summary statistics.
    """

    entity: str = ""
    total_number_of_entities_processed: int = 0
    number_of_unique_activity_units: int = 0
    number_of_unique_activity_units_enabled_by_gene_product_association: int = 0
    number_of_unique_activity_units_enabled_by_protein_complex_association: int = 0

    number_of_genes: int = 0
    number_of_unique_gene_product_enablers: int = 0
    number_of_unique_protein_complex_terms: int = 0
    number_of_unique_member_protein_complex_genes: int = 0
    number_of_unique_gene_product_and_protein_complex_gene_enablers: int = 0

    average_number_of_unique_gene_product_enablers_for_entity: float = 0.0

    number_of_causal_relations: int = 0
    number_of_unique_references: int = 0
    number_of_unique_pmid: int = 0
    number_of_inputs: int = 0
    number_of_chemical_inputs: int = 0
    number_of_other_inputs: int = 0
    number_of_unique_chemical_inputs: int = 0
    number_of_unique_other_inputs: int = 0
    number_of_outputs: int = 0
    number_of_chemical_outputs: int = 0
    number_of_other_outputs: int = 0
    number_of_unique_chemical_outputs: int = 0
    number_of_unique_other_outputs: int = 0
    number_of_go_terms: int = 0
    number_of_unique_go_terms: int = 0
    number_of_unique_inferred_relations: int = 0
    list_inferred_relations: List[List[str]] | None = []

    set_activities: Set[str] = set()
    set_enabled_by_gene_product: Set[str] = set()
    set_references: Set[str] = set()
    set_activity_unit_gene_product_enablers: Set[str] = set()
    set_activity_unit_protein_complex_enablers: Set[str] = set()
    set_protein_complex_in_activity_term: Set[str] = set()
    set_protein_complex_genes: Set[str] = set()
    list_model_details: List[ModelDetails] | None = []
    list_has_input_term: (
        List[str] | None
    ) = []  # Covers both has_input and has primary_input
    list_has_output_term: (
        List[str] | None
    ) = []  # Covers both has_output and has primary_output
    list_go_terms: List[str] | None = []


#: Return type for :func:`process_gocam_model_file`.
#: A 3-tuple of (result_type, query_index_or_none, reason_or_none).
ProcessingResult: TypeAlias = (
    tuple[Literal[ResultType.SUCCESS], QueryIndex, None]
    | tuple[Literal[ResultType.FILTERED], QueryIndex, FilterReason]
    | tuple[Literal[ResultType.ERROR], None, ErrorReason]
)


def setup_logger(verbose: int) -> None:
    """Set up the logger with the specified verbosity level.

    Args:
        verbose (int): Verbosity level (0: WARNING, 1: INFO, 2: DEBUG)
    """
    if verbose == 0:
        level = logging.WARNING
    elif verbose == 1:
        level = logging.INFO
    else:
        level = logging.DEBUG
    logging.basicConfig(level=level, format="%(message)s", handlers=[RichHandler()])


def _count_pmids(references: Collection[str]) -> int:
    """Count the number of PMID references in a set of references."""
    return count_type_in_collection(references, "pmid")


def count_chebis(terms: Collection[str]) -> int:
    """Count the number of CHEBI terms"""
    return count_type_in_collection(terms, "chebi")


def count_type_in_collection(string_col: Collection[str], prefix: str) -> int:
    """Count strings in a collection whose lowercase form starts with the given prefix.

    Args:
        string_col: Collection of strings to search.
        prefix: Case-insensitive prefix to match against.

    Returns:
        Number of matching strings.
    """
    return sum(1 for ent in string_col if ent.lower().startswith(prefix))


def _build_obsolete_ids(gocam_model: Model) -> set[str]:
    """Build a set of object IDs marked as obsolete in the model."""
    obsolete_ids: set[str] = set()
    for obj in gocam_model.objects or []:
        if getattr(obj, "obsolete", None) is True:
            obsolete_ids.add(obj.id)
    return obsolete_ids


def _is_obsolete(obj) -> bool:
    """Check if an object has obsolete=True."""
    return getattr(obj, "obsolete", None) is True


def compute_molecule_and_term_counts(stats: GocamStats | AggregateInfo) -> None:
    """Compute input/output/GO-term counts from the raw term lists on a stats object."""
    # Inputs
    input_terms = stats.list_has_input_term or []
    stats.number_of_inputs = len(input_terms)
    stats.number_of_chemical_inputs = count_chebis(input_terms)
    unique_inputs = set(input_terms)
    stats.number_of_unique_chemical_inputs = count_chebis(unique_inputs)
    stats.number_of_other_inputs = (
        stats.number_of_inputs - stats.number_of_chemical_inputs
    )
    stats.number_of_unique_other_inputs = (
        len(unique_inputs) - stats.number_of_unique_chemical_inputs
    )

    # Outputs
    output_terms = stats.list_has_output_term or []
    stats.number_of_outputs = len(output_terms)
    stats.number_of_chemical_outputs = count_chebis(output_terms)
    unique_outputs = set(output_terms)
    stats.number_of_unique_chemical_outputs = count_chebis(unique_outputs)
    stats.number_of_other_outputs = (
        stats.number_of_outputs - stats.number_of_chemical_outputs
    )
    stats.number_of_unique_other_outputs = (
        len(unique_outputs) - stats.number_of_unique_chemical_outputs
    )

    # GO terms
    go_terms = stats.list_go_terms or []
    stats.number_of_go_terms = len(go_terms)
    stats.number_of_unique_go_terms = len(set(go_terms))


def _update_entity_gene_stats(
    entity_info: GocamStats,
    activity,
    gocam_model_id: str,
    obsolete_ids: set[str],
) -> None:
    """Update gene product and protein complex enabler stats for an entity.

    Args:
        entity_info: The GocamStats object for the contributor or provider.
        activity: An Activity from the GO-CAM model.
        gocam_model_id: Identifier of the parent model.
        obsolete_ids: Set of object IDs marked as obsolete.
    """
    if isinstance(activity.enabled_by, EnabledByGeneProductAssociation):
        entity_info.set_activity_unit_gene_product_enablers.add(activity.id)
        if activity.enabled_by.term and activity.enabled_by.term not in obsolete_ids:
            if entity_info.list_enabled_by_gene_product is not None:
                entity_info.list_enabled_by_gene_product.append(
                    activity.enabled_by.term
                )
            entity_info.set_enabled_by_gene_product.add(activity.enabled_by.term)
            entity_info.number_of_genes += 1
    if isinstance(activity.enabled_by, EnabledByProteinComplexAssociation):
        entity_info.set_activity_unit_protein_complex_enablers.add(activity.id)
        if activity.enabled_by.term and activity.enabled_by.term not in obsolete_ids:
            entity_info.set_protein_complex_in_activity_term.add(
                activity.enabled_by.term
            )
        if activity.enabled_by.members:
            for member in activity.enabled_by.members:
                if member.term and member.term not in obsolete_ids:
                    entity_info.set_protein_complex_genes.add(member.term)
                    entity_info.number_of_genes += 1
    entity_info.set_activities.add(activity.id)
    entity_info.set_models.add(gocam_model_id)


def _is_association_obsolete(assoc: Association, obsolete_ids: set[str]) -> bool:
    """Check if an association references an obsolete term or molecule."""
    term = getattr(assoc, "term", None) or getattr(assoc, "molecule", None)
    return isinstance(term, str) and term in obsolete_ids


def _iter_activity_associations(
    activity,
    obsolete_ids: set[str],
) -> list[Association]:
    """Collect all non-obsolete Association objects from an Activity.

    Walks enabled_by, molecular_function, part_of, occurs_in, inputs,
    outputs, and causal_associations, skipping any that reference obsolete terms.

    Args:
        activity: An Activity from the GO-CAM model.
        obsolete_ids: Set of object IDs marked as obsolete.

    Returns:
        List of non-obsolete Association objects.
    """
    associations: list[Association] = []
    if activity.enabled_by and not _is_association_obsolete(
        activity.enabled_by, obsolete_ids
    ):
        associations.append(activity.enabled_by)
        if isinstance(activity.enabled_by, EnabledByProteinComplexAssociation):
            if activity.enabled_by.members:
                for member in activity.enabled_by.members:
                    if not _is_association_obsolete(member, obsolete_ids):
                        associations.append(member)
    if activity.molecular_function and not _is_association_obsolete(
        activity.molecular_function, obsolete_ids
    ):
        associations.append(activity.molecular_function)
    if activity.part_of and not _is_association_obsolete(
        activity.part_of, obsolete_ids
    ):
        associations.append(activity.part_of)
        if activity.part_of.happens_during and not _is_association_obsolete(
            activity.part_of.happens_during, obsolete_ids
        ):
            associations.append(activity.part_of.happens_during)
        if activity.part_of.part_of and not _is_association_obsolete(
            activity.part_of.part_of, obsolete_ids
        ):
            associations.append(activity.part_of.part_of)
    if activity.occurs_in and not _is_association_obsolete(
        activity.occurs_in, obsolete_ids
    ):
        associations.append(activity.occurs_in)
        if activity.occurs_in.part_of and not _is_association_obsolete(
            activity.occurs_in.part_of, obsolete_ids
        ):
            associations.append(activity.occurs_in.part_of)
    for mol_assoc in activity.molecular_associations or []:
        if not _is_association_obsolete(mol_assoc, obsolete_ids):
            associations.append(mol_assoc)
    for causal_association in activity.causal_associations or []:
        if not _is_association_obsolete(causal_association, obsolete_ids):
            associations.append(causal_association)
    return associations


def _collect_references(
    gocam_model: Model,
    stats_by_model: GocamStats,
    model_aggregate: AggregateInfo,
    contributor_lookup: Dict[str, GocamStats],
    provider_lookup: Dict[str, GocamStats],
    obsolete_ids: set[str],
) -> None:
    """Collect reference objects from all evidence items in the model.

    Iterates over every association across all activities, extracts references
    from evidence items, and attributes them to contributors and providers
    via the evidence-level provenances.

    Populates:
        stats_by_model.set_references: all references for this model
        model_aggregate.set_references: accumulated references across models
        contributor_lookup[contributor].set_references: references by contributor
        provider_lookup[provider].set_references: references by provider
    """
    for activity in gocam_model.activities or []:
        for association in _iter_activity_associations(activity, obsolete_ids):
            if not association.evidence:
                continue
            for evidence in association.evidence:
                if not evidence.reference:
                    continue
                ref = evidence.reference
                stats_by_model.set_references.add(ref)
                model_aggregate.set_references.add(ref)
                if evidence.provenances:
                    for provenance in evidence.provenances:
                        if provenance.contributor:
                            for contributor in provenance.contributor:
                                contributor_info = contributor_lookup.setdefault(
                                    contributor, GocamStats()
                                )
                                contributor_info.set_references.add(ref)
                        if provenance.provided_by:
                            for provider in provenance.provided_by:
                                provider_info = provider_lookup.setdefault(
                                    provider, GocamStats()
                                )
                                provider_info.set_references.add(ref)


def _collect_molecule_terms(
    activity,
    stats_by_model: GocamStats,
    model_aggregate: AggregateInfo,
    contributor_lookup: Dict[str, GocamStats],
    provider_lookup: Dict[str, GocamStats],
    obsolete_ids: set[str],
    molecule_lookup: Dict[str, str] | None = None,
) -> None:
    """Collect molecule input/output terms from MoleculeAssociation objects.

    Iterates over molecular_associations for the given activity, classifying
    each by its predicate into input or output, and attributes terms to
    contributors and providers via association-level provenances.

    Populates:
        stats_by_model.list_has_input_term / list_has_output_term
        model_aggregate.list_has_input_term / list_has_output_term
        contributor_lookup[contributor].list_has_input_term / list_has_output_term
        provider_lookup[provider].list_has_input_term / list_has_output_term
    """
    for ma in activity.molecular_associations or []:
        molecule_id = ma.molecule
        if not molecule_id or molecule_id in obsolete_ids:
            continue
        # Resolve molecule node ID to its actual term (e.g. CHEBI:58199)
        molecule = (
            molecule_lookup.get(molecule_id, molecule_id)
            if molecule_lookup
            else molecule_id
        )
        if not molecule or molecule in obsolete_ids:
            continue

        # Determine whether this is an input or output association
        if ma.predicate in INPUT_PREDICATES:
            target_attr = "list_has_input_term"
        elif ma.predicate in OUTPUT_PREDICATES:
            target_attr = "list_has_output_term"
        else:
            continue

        for target in (stats_by_model, model_aggregate):
            term_list = getattr(target, target_attr)
            if term_list is not None:
                term_list.append(molecule)

        if ma.provenances:
            for provenance in ma.provenances:
                if provenance.contributor:
                    for contributor in provenance.contributor:
                        contributor_info = contributor_lookup.setdefault(
                            contributor, GocamStats()
                        )
                        term_list = getattr(contributor_info, target_attr)
                        if term_list is not None:
                            term_list.append(molecule)
                if provenance.provided_by:
                    for provider in provenance.provided_by:
                        provider_info = provider_lookup.setdefault(
                            provider, GocamStats()
                        )
                        term_list = getattr(provider_info, target_attr)
                        if term_list is not None:
                            term_list.append(molecule)


def _collect_terms(
    activity,
    stats_by_model: GocamStats,
    model_aggregate: AggregateInfo,
    contributor_lookup: Dict[str, GocamStats],
    provider_lookup: Dict[str, GocamStats],
    obsolete_ids: set[str],
) -> None:
    """Collect GO terms from all objects within an activity.

    Recursively walks the entire activity object tree.  For every object
    that has a ``term`` attribute whose value starts with "GO:"
    (case-insensitive), the term is appended to list_go_terms on the
    model-level, aggregate, contributor, and provider stats objects.
    Contributor/provider attribution uses the ``provenances`` on the same
    object that carries the term.

    Populates:
        stats_by_model.list_go_terms
        model_aggregate.list_go_terms
        contributor_lookup[contributor].list_go_terms
        provider_lookup[provider].list_go_terms
    """

    def _walk(obj) -> None:
        if obj is None:
            return
        if not isinstance(obj, ConfiguredBaseModel):
            return

        # Skip obsolete objects and their subtrees
        if _is_obsolete(obj):
            return

        # Check whether this object carries a GO term
        term = getattr(obj, "term", None)
        if isinstance(term, str) and term in obsolete_ids:
            return
        if isinstance(term, str) and term.upper().startswith("GO:"):
            if stats_by_model.list_go_terms is not None:
                stats_by_model.list_go_terms.append(term)
            if model_aggregate.list_go_terms is not None:
                model_aggregate.list_go_terms.append(term)
            provenances = getattr(obj, "provenances", None)
            if provenances:
                for provenance in provenances:
                    if provenance.contributor:
                        for contributor in provenance.contributor:
                            contributor_info = contributor_lookup.setdefault(
                                contributor, GocamStats()
                            )
                            if contributor_info.list_go_terms is not None:
                                contributor_info.list_go_terms.append(term)
                    if provenance.provided_by:
                        for provider in provenance.provided_by:
                            provider_info = provider_lookup.setdefault(
                                provider, GocamStats()
                            )
                            if provider_info.list_go_terms is not None:
                                provider_info.list_go_terms.append(term)

        # Recurse into all fields of this Pydantic model
        for field_name in type(obj).model_fields:
            value = getattr(obj, field_name, None)
            if value is None:
                continue
            if isinstance(value, list):
                for item in value:
                    _walk(item)
            else:
                _walk(value)

    _walk(activity)


def _get_chemical_terms(
    molecular_associations: list | None,
    predicates: set[str],
    obsolete_ids: set[str],
    molecule_lookup: Dict[str, str] | None = None,
) -> set[str]:
    """Extract CHEBI chemical terms from molecular associations matching given predicates.

    Filters molecular_associations by predicate, then collects molecule values
    that are CHEBI terms and not obsolete.

    Args:
        molecular_associations: List of MoleculeAssociation objects from an activity.
        predicates: Set of RO relation URIs to match (e.g. INPUT_PREDICATES or OUTPUT_PREDICATES).
        obsolete_ids: Set of object IDs marked as obsolete.

    Returns:
        A set of CHEBI molecule strings.
    """
    terms: set[str] = set()
    for ma in molecular_associations or []:
        if ma.predicate not in predicates:
            continue
        molecule_id = ma.molecule
        if not molecule_id or molecule_id in obsolete_ids:
            continue
        # Resolve molecule node ID to its actual term (e.g. CHEBI:58199)
        molecule = (
            molecule_lookup.get(molecule_id, molecule_id)
            if molecule_lookup
            else molecule_id
        )
        if (
            molecule
            and molecule not in obsolete_ids
            and molecule.lower().startswith("chebi")
        ):
            terms.add(molecule)
    return terms


def _compute_inferred_relations(
    activities: list,
    obsolete_ids: set[str],
    molecule_lookup: Dict[str, str] | None = None,
) -> list[list[str]]:
    """Compute inferred relations between activities based on shared chemicals.

    An inferred relation exists when activity A produces chemical outputs that
    are consumed as inputs by activity B, and activity B also produces chemical
    outputs. This captures implicit metabolic pathway connections where one
    activity's product feeds into another activity's transformation.

    Args:
        activities: List of Activity objects from a GO-CAM model.
        obsolete_ids: Set of object IDs marked as obsolete.

    Returns:
        A deduplicated list of [activity_A_id, activity_B_id] pairs representing
        inferred relations.
    """
    # Build per-activity chemical input and output sets
    activity_chem_outputs: dict[str, set[str]] = {}
    activity_chem_inputs: dict[str, set[str]] = {}

    for activity in activities:
        outputs = _get_chemical_terms(
            activity.molecular_associations,
            OUTPUT_PREDICATES,
            obsolete_ids,
            molecule_lookup,
        )
        if outputs:
            activity_chem_outputs[activity.id] = outputs

        inputs = _get_chemical_terms(
            activity.molecular_associations,
            INPUT_PREDICATES,
            obsolete_ids,
            molecule_lookup,
        )
        if inputs:
            activity_chem_inputs[activity.id] = inputs

    # Find pairs where A's outputs overlap with B's inputs and B also has outputs
    relations: list[list[str]] = []
    seen: set[tuple[str, str]] = set()

    for a_id, a_outputs in activity_chem_outputs.items():
        for b_id, b_inputs in activity_chem_inputs.items():
            if a_id == b_id:
                continue
            if a_outputs & b_inputs and b_id in activity_chem_outputs:
                pair = (a_id, b_id)
                if pair not in seen:
                    seen.add(pair)
                    relations.append([a_id, b_id])

    return relations


def process_gocam_model_file(
    json_file: Path,
    output_dir: Path | None,
    model_aggregate: AggregateInfo,
    contributor_lookup: Dict[str, GocamStats],
    provider_lookup: Dict[str, GocamStats],
) -> ProcessingResult:
    """Process a single GO-CAM model JSON file and compute statistics.

    Reads the model, indexes it, collects per-activity counts (enablers,
    causal relations, molecule terms, GO terms, inferred relations, and
    references), and writes per-model statistics files when an output
    directory is provided.

    Args:
        json_file: Path to the GO-CAM model JSON file.
        output_dir: Directory to save the GO-CAM model statistics files.
            If None (dry run), no files are written.
        model_aggregate: Accumulator for cross-model aggregate statistics.
        contributor_lookup: Mapping of contributor URI to per-contributor stats.
        provider_lookup: Mapping of provider URI to per-provider stats.

    Returns:
        A ProcessingResult 3-tuple of (result_type, query_index_or_none,
        reason_or_none).
    """
    indexer = Indexer()
    # Read GO-CAM JSON file
    try:
        with open(json_file, "r") as f:
            # json_content = f.read()
            # data = json.loads(json_content)
            # model_json = json.dumps(data)
            # gocam_model = Model.model_validate_json(model_json)
            gocam_model = Model.model_validate_json(f.read())
        logger.debug(f"Successfully read GO-CAM model from {json_file}")
    except Exception as e:
        logger.error(f"Error reading file {json_file}", exc_info=e)
        return ResultType.ERROR, None, ErrorReason.READ_ERROR

    # Populate the model with indexing information
    indexer.index_model(gocam_model)

    # Build set of obsolete object IDs to filter from statistics
    obsolete_ids = _build_obsolete_ids(gocam_model)

    # Build molecule node ID → term lookup (e.g. "gomodel:.../id" → "CHEBI:58199")
    molecule_lookup: Dict[str, str] = {}
    for mol_node in gocam_model.molecules or []:
        if mol_node.id and mol_node.term:
            molecule_lookup[mol_node.id] = mol_node.term

    # Get model statistics information
    calculated_aggregate_values_by_model = gocam_model.query_index

    # Detailed model statistics information
    stats_by_model = GocamStats()
    model_details = ModelDetails()
    model_details.file_name = json_file.name
    model_details.model_id = gocam_model.id
    model_details.model_name = gocam_model.title
    if model_aggregate.list_model_details is not None:
        model_aggregate.list_model_details.append(model_details)
    stats_by_model.number_of_models = 1

    if gocam_model.activities:
        for activity in gocam_model.activities:
            stats_by_model.set_activities.add(activity.id)
            model_aggregate.set_activities.add(activity.id)
            model_details.set_activities.add(activity.id)
            if isinstance(activity.enabled_by, EnabledByGeneProductAssociation):
                stats_by_model.set_activity_unit_gene_product_enablers.add(activity.id)
                model_aggregate.set_activity_unit_gene_product_enablers.add(activity.id)
                if (
                    activity.enabled_by.term
                    and activity.enabled_by.term not in obsolete_ids
                ):
                    if stats_by_model.list_enabled_by_gene_product is not None:
                        stats_by_model.list_enabled_by_gene_product.append(
                            activity.enabled_by.term
                        )
                    stats_by_model.set_enabled_by_gene_product.add(
                        activity.enabled_by.term
                    )
                    stats_by_model.number_of_genes += 1
                    model_aggregate.set_enabled_by_gene_product.add(
                        activity.enabled_by.term
                    )

            elif isinstance(activity.enabled_by, EnabledByProteinComplexAssociation):
                stats_by_model.set_activity_unit_protein_complex_enablers.add(
                    activity.id
                )
                model_aggregate.set_activity_unit_protein_complex_enablers.add(
                    activity.id
                )
                if (
                    activity.enabled_by.term
                    and activity.enabled_by.term not in obsolete_ids
                ):
                    stats_by_model.set_protein_complex_in_activity_term.add(
                        activity.enabled_by.term
                    )
                    model_aggregate.set_protein_complex_in_activity_term.add(
                        activity.enabled_by.term
                    )
                if activity.enabled_by.members:
                    for member in activity.enabled_by.members:
                        if member.term and member.term not in obsolete_ids:
                            stats_by_model.set_protein_complex_genes.add(member.term)
                            model_aggregate.set_protein_complex_genes.add(member.term)
                            stats_by_model.number_of_genes += 1

            if (
                activity.enabled_by is not None
                and activity.enabled_by.provenances is not None
            ):
                for provenance in activity.enabled_by.provenances:
                    if provenance.contributor:
                        for contributor in provenance.contributor:
                            contributor_info = contributor_lookup.setdefault(
                                contributor, GocamStats()
                            )
                            _update_entity_gene_stats(
                                contributor_info,
                                activity,
                                gocam_model.id,
                                obsolete_ids,
                            )
                    if provenance.provided_by:
                        for provider in provenance.provided_by:
                            provider_info = provider_lookup.setdefault(
                                provider, GocamStats()
                            )
                            _update_entity_gene_stats(
                                provider_info,
                                activity,
                                gocam_model.id,
                                obsolete_ids,
                            )

            # Add information about causal relations
            if activity.causal_associations:
                for causal_association in activity.causal_associations:
                    if stats_by_model.list_causal_relations is not None:
                        stats_by_model.list_causal_relations.append(activity.id)
                    stats_by_model.set_causal_relations.add(activity.id)
                    if causal_association.provenances:
                        for provenance in causal_association.provenances:
                            if provenance.contributor:
                                for contributor in provenance.contributor:
                                    contributor_info = contributor_lookup.setdefault(
                                        contributor, GocamStats()
                                    )
                                    contributor_info.set_causal_relations.add(
                                        activity.id
                                    )
                                    contributor_info.number_of_causal_relations += 1
                            if provenance.provided_by:
                                for provider in provenance.provided_by:
                                    provider_info = provider_lookup.setdefault(
                                        provider, GocamStats()
                                    )
                                    provider_info.set_causal_relations.add(activity.id)
                                    provider_info.number_of_causal_relations += 1

            # Handle has input and has output
            _collect_molecule_terms(
                activity,
                stats_by_model,
                model_aggregate,
                contributor_lookup,
                provider_lookup,
                obsolete_ids,
                molecule_lookup,
            )

            # Collect GO terms from all associations
            _collect_terms(
                activity,
                stats_by_model,
                model_aggregate,
                contributor_lookup,
                provider_lookup,
                obsolete_ids,
            )

    # Collect references from all evidence items across all associations
    _collect_references(
        gocam_model,
        stats_by_model,
        model_aggregate,
        contributor_lookup,
        provider_lookup,
        obsolete_ids,
    )

    # Compute inferred relations for the model
    if gocam_model.activities:
        inferred = _compute_inferred_relations(
            gocam_model.activities, obsolete_ids, molecule_lookup
        )
        stats_by_model.list_inferred_relations = inferred
        stats_by_model.number_of_unique_inferred_relations = len(inferred)
        if model_aggregate.list_inferred_relations is not None:
            model_aggregate.list_inferred_relations.extend(inferred)

    # Set fields, counts and sort data for current model
    stats_by_model.number_of_causal_relations = len(
        stats_by_model.list_causal_relations or []
    )
    stats_by_model.number_of_unique_causal_relations = len(
        stats_by_model.set_causal_relations
    )
    stats_by_model.number_of_activity_units = (
        calculated_aggregate_values_by_model.number_of_activities
    )
    stats_by_model.number_of_unique_references = len(stats_by_model.set_references)
    stats_by_model.number_of_unique_pmid = _count_pmids(stats_by_model.set_references)

    stats_by_model.number_of_activity_units_enabled_by_gene_product = len(
        stats_by_model.set_activity_unit_gene_product_enablers
    )
    stats_by_model.number_of_activity_units_enabled_by_protein_complex = len(
        stats_by_model.set_activity_unit_protein_complex_enablers
    )

    stats_by_model.number_of_unique_gene_product_enablers = len(
        stats_by_model.set_enabled_by_gene_product
    )
    stats_by_model.number_of_unique_protein_complex_genes = len(
        stats_by_model.set_protein_complex_genes
    )
    stats_by_model.number_of_unique_gene_product_and_protein_complex_gene_enablers = (
        len(
            stats_by_model.set_enabled_by_gene_product.union(
                stats_by_model.set_protein_complex_genes
            )
        )
    )
    if stats_by_model.list_enabled_by_gene_product is not None:
        stats_by_model.list_enabled_by_gene_product.sort()
    compute_molecule_and_term_counts(stats_by_model)

    # Update aggregate model data
    model_aggregate.total_number_of_entities_processed += 1
    model_aggregate.number_of_genes = (
        model_aggregate.number_of_genes + stats_by_model.number_of_genes
    )
    model_aggregate.number_of_causal_relations = (
        model_aggregate.number_of_causal_relations
        + stats_by_model.number_of_causal_relations
    )

    # If dry run is enabled, skip writing the output file
    if output_dir is None:
        logger.info(
            f"Dry run enabled; skipping write for GO-CAM model {gocam_model.id}"
        )
        return ResultType.SUCCESS, calculated_aggregate_values_by_model, None

    # Write GO-CAM stats to output directory
    calculated_aggregate_values_by_model_subdir = "calculated_aggregate_values_by_model"
    calculated_aggregate_values_by_model_name = (
        "calculated_aggregate_values_by model_" + json_file.name
    )
    calculated_aggregate_values_by_model_name_file = (
        output_dir
        / calculated_aggregate_values_by_model_subdir
        / calculated_aggregate_values_by_model_name
    )
    calculated_aggregate_values_by_model_name_file.parent.mkdir(
        parents=True, exist_ok=True
    )

    stats_by_model_subdir = "stats_by_model"
    stats_by_model_name = "stats_by_model_" + json_file.name
    stats_by_model_output_file = (
        output_dir / stats_by_model_subdir / stats_by_model_name
    )
    stats_by_model_output_file.parent.mkdir(parents=True, exist_ok=True)

    try:
        model_stats_json = calculated_aggregate_values_by_model.model_dump_json(
            exclude_none=True
        )
        with open(calculated_aggregate_values_by_model_name_file, "w") as f:
            f.write(model_stats_json)
        logger.info(
            f"Successfully wrote GO-CAM all_stats to {calculated_aggregate_values_by_model_name_file}"
        )

        stats_by_model_json = stats_by_model.model_dump_json(exclude_none=True)
        with open(stats_by_model_output_file, "w") as f:
            f.write(stats_by_model_json)
        logger.info(
            f"Successfully wrote GO-CAM detailed_model_stats to {stats_by_model_output_file}"
        )

    except Exception as e:
        logger.error(f"An exception has occurred: {e}")
        return ResultType.ERROR, None, ErrorReason.WRITE_ERROR

    # If we reach here, the conversion and writing were successful
    return ResultType.SUCCESS, calculated_aggregate_values_by_model, None


def create_filename_from_url(url: str, replacement: str = "_") -> str:
    """Generate a filesystem-safe filename from a URL.

    Extracts the path component, sanitizes invalid characters, and
    truncates to 250 characters.

    Args:
        url: The URL to derive a filename from.
        replacement: Character used to replace invalid filename characters.

    Returns:
        A sanitized filename string.
    """
    # 1. Parse the URL to get the path component
    parsed_url = urlparse(url)
    path = parsed_url.path

    # 2. Extract the base filename from the path
    # Use os.path.basename for handling different path separators, though urlparse usually handles this
    filename = os.path.basename(path)

    # If the URL ends with a '/', filename might be empty, use the hostname as a fallback
    if not filename:
        filename = parsed_url.hostname or "downloaded_file"

    # 3. Sanitize the filename for cross-platform compatibility
    # Invalid characters are typically /\\:*?"<>| and null character
    safe_filename = re.sub(r'[\\/:*?"<>| \n\t]', replacement, filename)

    # 4. Truncate if necessary (some file systems have length limits, e.g., 255 characters)
    max_length = 250
    if len(safe_filename) > max_length:
        safe_filename = safe_filename[:max_length]

    # Ensure filename is not empty or invalid (e.g., "." or "..")
    if not safe_filename or safe_filename in (".", ".."):
        safe_filename = "safe_download"

    return safe_filename


def output_entity_results(
    output_dir: Path | None,
    entity_label: str,
    entity_lookup: Dict[str, GocamStats],
    entity_sub_dir: str,
    entity_agg_file_name: str,
) -> None:
    """Compute final counts and write per-entity and aggregate statistics.

    Iterates over all entities (contributors or providers), finalizes their
    numeric counts from working sets, accumulates into an entity-level
    aggregate, and writes JSON output files.

    Args:
        output_dir: Directory to write JSON files. If None, no files are written.
        entity_label: Human-readable label for this entity type (e.g. "Curator", "Group").
        entity_lookup: Mapping of entity URI to per-entity GocamStats.
        entity_sub_dir: Subdirectory name for per-entity output files.
        entity_agg_file_name: Filename for the aggregate statistics JSON.
    """
    # Inferred relations are only meaningful at the model level, not per-entity
    _exclude_inferred = {
        "number_of_unique_inferred_relations": True,
        "list_inferred_relations": True,
    }

    entity_agg = AggregateInfo()
    entity_agg.entity = entity_label
    entity_agg.total_number_of_entities_processed = len(entity_lookup)

    for entity, details in entity_lookup.items():
        logger.debug(f"Processing contributor: {entity}")
        # Set numbers
        details.uri = entity
        details.number_of_models = len(details.set_models)
        details.number_of_unique_gene_product_enablers = len(
            details.set_enabled_by_gene_product
        )
        details.number_of_unique_protein_complex_genes = len(
            details.set_protein_complex_genes
        )
        details.number_of_unique_gene_product_and_protein_complex_gene_enablers = len(
            details.set_enabled_by_gene_product.union(details.set_protein_complex_genes)
        )
        details.number_of_unique_references = len(details.set_references)
        details.number_of_activity_units = len(details.set_activities)
        details.number_of_activity_units_enabled_by_gene_product = len(
            details.set_activity_unit_gene_product_enablers
        )
        details.number_of_activity_units_enabled_by_protein_complex = len(
            details.set_activity_unit_protein_complex_enablers
        )
        details.number_of_unique_causal_relations = len(details.set_causal_relations)
        compute_molecule_and_term_counts(details)

        # Sort lists
        if details.list_enabled_by_gene_product is not None:
            details.list_enabled_by_gene_product.sort()
        details.number_of_unique_pmid = _count_pmids(details.set_references)

        entity_agg.number_of_genes = (
            entity_agg.number_of_genes + details.number_of_genes
        )
        entity_agg.number_of_causal_relations = (
            entity_agg.number_of_causal_relations + details.number_of_causal_relations
        )
        entity_agg.set_activities.update(details.set_activities)
        entity_agg.set_enabled_by_gene_product.update(
            details.set_enabled_by_gene_product
        )
        entity_agg.set_references.update(details.set_references)
        entity_agg.set_activity_unit_gene_product_enablers.update(
            details.set_activity_unit_gene_product_enablers
        )
        entity_agg.set_activity_unit_protein_complex_enablers.update(
            details.set_activity_unit_protein_complex_enablers
        )
        entity_agg.set_protein_complex_in_activity_term.update(
            details.set_protein_complex_in_activity_term
        )
        entity_agg.set_protein_complex_genes.update(details.set_protein_complex_genes)
        if (
            entity_agg.list_has_input_term is not None
            and details.list_has_input_term is not None
        ):
            entity_agg.list_has_input_term.extend(details.list_has_input_term)
        if (
            entity_agg.list_has_output_term is not None
            and details.list_has_output_term is not None
        ):
            entity_agg.list_has_output_term.extend(details.list_has_output_term)
        if entity_agg.list_go_terms is not None and details.list_go_terms is not None:
            entity_agg.list_go_terms.extend(details.list_go_terms)
        compute_molecule_and_term_counts(entity_agg)

        if output_dir is not None:
            stats_by_entity_subdir = entity_sub_dir
            json_file_name = create_filename_from_url(entity + ".json")
            stats_by_curator_name = entity_sub_dir + "_" + json_file_name
            stats_by_curator_output_file = (
                output_dir / stats_by_entity_subdir / stats_by_curator_name
            )
            stats_by_curator_output_file.parent.mkdir(parents=True, exist_ok=True)

            try:
                contributor_model_json = details.model_dump_json(
                    exclude_none=True, exclude=_exclude_inferred
                )
                with open(stats_by_curator_output_file, "w") as f:
                    f.write(contributor_model_json)
                logger.info(
                    f"Successfully wrote contributor model to {stats_by_curator_output_file}"
                )
            except Exception as e:
                logger.error(
                    f"Error writing contributor model file {stats_by_curator_output_file}",
                    exc_info=e,
                )
                return ResultType.ERROR, ErrorReason.WRITE_ERROR

    entity_agg.number_of_unique_activity_units = len(entity_agg.set_activities)
    entity_agg.number_of_unique_gene_product_enablers = len(
        entity_agg.set_enabled_by_gene_product
    )
    entity_agg.number_of_unique_references = len(entity_agg.set_references)
    entity_agg.number_of_unique_pmid = _count_pmids(entity_agg.set_references)
    entity_agg.number_of_unique_activity_units_enabled_by_protein_complex_association = len(
        entity_agg.set_activity_unit_protein_complex_enablers
    )
    entity_agg.number_of_unique_activity_units_enabled_by_gene_product_association = (
        len(entity_agg.set_activity_unit_gene_product_enablers)
    )
    entity_agg.number_of_unique_protein_complex_terms = len(
        entity_agg.set_protein_complex_in_activity_term
    )
    entity_agg.number_of_unique_member_protein_complex_genes = len(
        entity_agg.set_protein_complex_genes
    )
    entity_agg.number_of_unique_gene_product_and_protein_complex_gene_enablers = len(
        entity_agg.set_enabled_by_gene_product.union(
            entity_agg.set_protein_complex_genes
        )
    )

    if entity_agg.total_number_of_entities_processed != 0:
        entity_agg.average_number_of_unique_gene_product_enablers_for_entity = round(
            entity_agg.number_of_unique_gene_product_enablers
            / entity_agg.total_number_of_entities_processed,
            2,
        )

    if output_dir is not None:
        aggregate_stats_name = entity_agg_file_name
        aggregate_stats_name_output_file = output_dir / aggregate_stats_name
        aggregate_stats_name_output_file.parent.mkdir(parents=True, exist_ok=True)

        try:
            aggregate_json = entity_agg.model_dump_json(
                exclude_none=True, exclude=_exclude_inferred
            )
            with open(aggregate_stats_name_output_file, "w") as f:
                f.write(aggregate_json)
            logger.info(
                f"Successfully wrote entity aggregate stats to {aggregate_stats_name_output_file}"
            )
        except Exception as e:
            logger.error(
                f"Error writing entity aggregate stats {aggregate_stats_name_output_file}",
                exc_info=e,
            )
            return ResultType.ERROR, ErrorReason.WRITE_ERROR


def output_summary(
    results: list[tuple[Path, ProcessingResult]],
    output_dir: Path | None,
    model_aggregate: AggregateInfo,
    contributor_lookup: Dict[str, GocamStats],
    provider_lookup: Dict[str, GocamStats],
) -> None:
    """Finalize aggregate statistics and write all summary output files.

    Computes final counts on the model-level aggregate, writes the
    aggregate model stats file, delegates contributor and provider
    aggregate output, and prints a success/failure summary tree.

    Args:
        results: List of (json_file, ProcessingResult) tuples from model processing.
        output_dir: Directory to write JSON files. If None, no files are written.
        model_aggregate: Accumulated cross-model aggregate statistics.
        contributor_lookup: Mapping of contributor URI to per-contributor stats.
        provider_lookup: Mapping of provider URI to per-provider stats.
    """
    model_aggregate.entity = "Model"
    model_aggregate.number_of_unique_gene_product_enablers = len(
        model_aggregate.set_enabled_by_gene_product
    )
    model_aggregate.number_of_unique_references = len(model_aggregate.set_references)
    model_aggregate.number_of_unique_activity_units = len(
        model_aggregate.set_activities
    )
    model_aggregate.number_of_unique_activity_units_enabled_by_protein_complex_association = len(
        model_aggregate.set_activity_unit_protein_complex_enablers
    )
    model_aggregate.number_of_unique_activity_units_enabled_by_gene_product_association = len(
        model_aggregate.set_activity_unit_gene_product_enablers
    )
    model_aggregate.number_of_unique_protein_complex_terms = len(
        model_aggregate.set_protein_complex_in_activity_term
    )
    model_aggregate.number_of_unique_member_protein_complex_genes = len(
        model_aggregate.set_protein_complex_genes
    )
    model_aggregate.number_of_unique_gene_product_and_protein_complex_gene_enablers = (
        len(
            model_aggregate.set_enabled_by_gene_product.union(
                model_aggregate.set_protein_complex_genes
            )
        )
    )
    compute_molecule_and_term_counts(model_aggregate)

    model_aggregate.number_of_unique_inferred_relations = len(
        model_aggregate.list_inferred_relations or []
    )

    if model_aggregate.total_number_of_entities_processed != 0:
        model_aggregate.average_number_of_unique_gene_product_enablers_for_entity = (
            round(
                model_aggregate.number_of_unique_gene_product_enablers
                / model_aggregate.total_number_of_entities_processed,
                2,
            )
        )

    model_aggregate.number_of_unique_pmid = _count_pmids(model_aggregate.set_references)

    if output_dir is not None:
        aggregate_stats_name = "aggregate_model_stats.json"
        aggregate_stats_name_output_file = output_dir / aggregate_stats_name
        aggregate_stats_name_output_file.parent.mkdir(parents=True, exist_ok=True)

        try:
            aggregate_json = model_aggregate.model_dump_json(exclude_none=True)
            with open(aggregate_stats_name_output_file, "w") as f:
                f.write(aggregate_json)
            logger.info(
                f"Successfully wrote model aggregate stats to {aggregate_stats_name_output_file}"
            )
        except Exception as e:
            logger.error(
                f"Error writing model aggregate stats {aggregate_stats_name_output_file}",
                exc_info=e,
            )
            return ResultType.ERROR, ErrorReason.WRITE_ERROR

    # Output information for contributor and provider
    output_entity_results(
        output_dir=output_dir,
        entity_label="Curator",
        entity_lookup=contributor_lookup,
        entity_sub_dir="stats_by_curator",
        entity_agg_file_name="aggregate_curator_stats.json",
    )
    output_entity_results(
        output_dir=output_dir,
        entity_label="Group",
        entity_lookup=provider_lookup,
        entity_sub_dir="stats_by_group",
        entity_agg_file_name="aggregate_group_stats.json",
    )

    total_count = len(results)
    success_count = sum(
        1 for _, (result, _, _) in results if result == ResultType.SUCCESS
    )
    failure_count = total_count - success_count

    tree = Tree(f"[bold]Processed {total_count} models[/bold]")
    if success_count > 0:
        tree.add(
            f"Successfully output stats for  [b]{success_count}[/b] models",
            style="green",
        )

    if failure_count > 0:
        tree.add(
            f"Did not output stats for  [b]{failure_count}[/b] models", style="red"
        )
    print(tree)


@app.command()
def main(
    input_dir: Annotated[
        Path,
        typer.Option(
            exists=True,
            file_okay=False,
            dir_okay=True,
            readable=True,
            help="Directory containing GO-CAM model json files.",
        ),
    ],
    output_dir: Annotated[
        Path | None,
        typer.Option(
            exists=True,
            file_okay=False,
            dir_okay=True,
            writable=True,
            help="Directory to save statistics information in json format. Required unless --dry-run is used.",
        ),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option(
            help="If set, the statistics will be calculated, but no files will be written.",
        ),
    ] = False,
    verbose: Annotated[
        int,
        typer.Option(
            "--verbose",
            "-v",
            count=True,
            help="Increase verbosity level. Can be used multiple times.",
        ),
    ] = 0,
    limit: Annotated[
        int,
        typer.Option(
            help="Limit the number of models to process. 0 means no limit.",
        ),
    ] = 0,
):
    """CLI entry point: process GO-CAM model JSON files and output statistics."""
    setup_logger(verbose)
    # Validate output directory if not a dry run
    if not dry_run and output_dir is None:
        raise typer.BadParameter(
            "Output directory must be specified unless --dry-run is used."
        )

    # Get list of JSON files in the input directory, excluding hidden files, and sort them
    # for consistent processing order. Apply limit if specified.
    json_files = sorted(
        f for f in input_dir.glob("*.json") if not f.name.startswith(".")
    )
    if limit > 0:
        json_files = json_files[:limit]

    # If no JSON files found, log a warning and exit
    if not json_files:
        logger.warning(f"No JSON files found in the specified directory: {input_dir}")
        raise typer.Exit(code=1)

    # Process each JSON file
    results: list[tuple[Path, ProcessingResult]] = []
    model_aggregate = AggregateInfo()
    contributor_lookup: Dict[str, GocamStats] = {}
    provider_lookup: Dict[str, GocamStats] = {}

    for json_file in track(
        json_files, description="Processing GO-CAM models and calculating statistics..."
    ):
        logger.debug(f"Processing file: {json_file}")
        result = process_gocam_model_file(
            json_file,
            output_dir=output_dir,
            model_aggregate=model_aggregate,
            contributor_lookup=contributor_lookup,
            provider_lookup=provider_lookup,
        )
        results.append((json_file, result))

    # Print result
    output_summary(
        results,
        output_dir=output_dir,
        model_aggregate=model_aggregate,
        contributor_lookup=contributor_lookup,
        provider_lookup=provider_lookup,
    )


if __name__ == "__main__":
    app()
