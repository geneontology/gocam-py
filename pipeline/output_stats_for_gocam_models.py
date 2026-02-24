"""
Output statistics information about GO-CAM models
"""

import logging
import os
import re
from collections import defaultdict
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

from rich import print
from rich.logging import RichHandler
from rich.progress import track
from rich.tree import Tree
import typer

from gocam.datamodel import (
    Association,
    ConfiguredBaseModel,
    EnabledByGeneProductAssociation,
    EnabledByProteinComplexAssociation,
    Model,
    MoleculeAssociation,
    MoleculeTermObject,
    QueryIndex,
)
from gocam.indexing.Indexer import Indexer

app = typer.Typer()

logger = logging.getLogger(__name__)

class ResultType(str, Enum):
    SUCCESS = "Success"
    FILTERED = "Filtered"
    ERROR = "Error"

class FilterReason(str, Enum):
    NO_MODEL = "No model"


class ErrorReason(str, Enum):
    READ_ERROR = "Read error"
    STATS_ERROR = "Error calculating statistics"
    WRITE_ERROR = "Write error"


class ModelStats(ConfiguredBaseModel):
    number_of_activitiy_units: int = 0
    number_of_activity_units_enabled_by_gene_product: int = 0
    number_of_activity_units_enabled_by_protein_complex: int = 0
    number_of_unique_gene_product_enablers: int = 0
    number_of_unique_protein_complex_genes: int = 0
    number_of_unique_gene_product_and_protein_complex_gene_enablers: int = 0
    number_of_genes: int = 0
    number_of_unique_references: int = 0
    number_of_unique_pmid: int = 0
    number_of_causal_relations: int = 0
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
    
    

    set_activities: Set[str] = set()  #= Field(default_factory=set)
    set_enabled_by_gene_product: Set[str] = set()
    set_protein_complex_genes: Set[str] = set()
    list_enabled_by_gene_product: List[str] | None = []  
    set_references: Set[str] = set()
    list_causal_relations: List[str] | None = []  
    set_gene_product_association_activity: Set[str] = set()
    set_protein_complex_association_activity: Set[str] = set()
    set_activity_unit_gene_product_enablers: Set[str] = set()
    set_activity_unit_protein_complex_enablers: Set[str] = set()
    set_protein_complex_in_activity_term: Set[str] = set()
    list_has_input_term:List[str] | None = [] # Covers both has_input and has primary_input
    set_has_input_term:List[str] | None = []
    list_has_output_term:List[str] | None = [] # Covers both has_output and has primary_output
    set_has_output_term:List[str] | None = []
    list_go_terms:List[str] | None = []


class ContributorProviderStats(ConfiguredBaseModel):
    number_of_models: int = 0
    number_of_activitiy_units: int = 0
    number_of_unique_activity_units_enabled_by_gene_product_association: int = 0
    number_of_unique_activity_units_enabled_by_protein_complex_association: int = 0
    number_of_unique_gene_product_enablers: int = 0
    number_of_unique_member_protein_complex_genes: int = 0
    number_of_unique_gene_product_and_protein_complex_gene_enablers: int = 0  

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

    set_models: Set[str] = set()
    set_activities: Set[str] = set()
    set_enabled_by_gene_product: Set[str] = set()
    list_enabled_by_gene_product: List[str]| None = []  
    set_references: Set[str] = set()
    set_causal_relations: Set[str] = set()
    set_activity_unit_gene_product_enablers: Set[str] = set()
    set_activity_unit_protein_complex_enablers: Set[str] = set()
    set_protein_complex_in_activity_term: Set[str] = set()
    set_protein_complex_genes: Set[str] = set()
    list_has_input_term:List[str] | None = [] # Covers both has_input and has primary_input
    set_has_input_term:List[str] | None = []
    list_has_output_term:List[str] | None = [] # Covers both has_output and has primary_output
    set_has_output_term:List[str] | None = []
    list_go_terms:List[str] | None = []

class ModelDetails(ConfiguredBaseModel):
    file_name: str = ""
    model_id: str = ""
    model_name: str = ""
    set_activities: Set[str] = set()

class AggregateInfo(ConfiguredBaseModel):
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

    set_activities: Set[str] = set()
    set_enabled_by_gene_product: Set[str] = set()
    set_references: Set[str] = set()
    set_activity_unit_gene_product_enablers: Set[str] = set()
    set_activity_unit_protein_complex_enablers: Set[str] = set()
    set_protein_complex_in_activity_term: Set[str] = set()
    set_protein_complex_genes: Set[str] = set()
    list_model_details: List[ModelDetails] | None = []
    list_has_input_term:List[str] | None = [] # Covers both has_input and has primary_input
    set_has_input_term:List[str] | None = []
    list_has_output_term:List[str] | None = [] # Covers both has_output and has primary_output
    set_has_output_term:List[str] | None = []
    list_go_terms:List[str] | None = []


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
    return sum(1 for ent in string_col if ent.lower().startswith(prefix))        


def _update_entity_evidence_stats(
    entity_info: ContributorProviderStats,
    activity,
    gocam_model_id: str,
) -> None:
    """Update contributor/provider stats from evidence data."""
    if isinstance(activity.enabled_by, EnabledByGeneProductAssociation):
        entity_info.set_activity_unit_gene_product_enablers.add(activity.id)
        if activity.enabled_by.term:
            entity_info.list_enabled_by_gene_product.append(activity.enabled_by.term)
            entity_info.set_enabled_by_gene_product.add(activity.enabled_by.term)
    if isinstance(activity.enabled_by, EnabledByProteinComplexAssociation):
        entity_info.set_activity_unit_protein_complex_enablers.add(activity.id)
        if activity.enabled_by.term:
            entity_info.set_protein_complex_in_activity_term.add(activity.enabled_by.term)
        if activity.enabled_by.members:
            for member in activity.enabled_by.members:
                if member.term:
                    entity_info.set_protein_complex_genes.add(member.term)
    entity_info.set_activities.add(activity.id)
    entity_info.set_models.add(gocam_model_id)


def _iter_activity_associations(activity) -> list[Association]:
    """Yield all Association objects from an Activity."""
    associations: list[Association] = []
    if activity.enabled_by:
        associations.append(activity.enabled_by)
        if isinstance(activity.enabled_by, EnabledByProteinComplexAssociation):
            if activity.enabled_by.members:
                associations.extend(activity.enabled_by.members)
    if activity.molecular_function:
        associations.append(activity.molecular_function)
    if activity.part_of:
        associations.append(activity.part_of)
        if activity.part_of.happens_during:
            associations.append(activity.part_of.happens_during)
        if activity.part_of.part_of:
            associations.append(activity.part_of.part_of)
    if activity.occurs_in:
        associations.append(activity.occurs_in)
        if activity.occurs_in.part_of:
            associations.append(activity.occurs_in.part_of)
    if activity.has_primary_input:
        associations.append(activity.has_primary_input)
    for has_input in activity.has_input or []:
        associations.append(has_input)
    if activity.has_primary_output:
        associations.append(activity.has_primary_output)
    for has_output in activity.has_output or []:
        associations.append(has_output)
    for causal_association in activity.causal_associations or []:
        associations.append(causal_association)
    return associations


def _collect_references(
    gocam_model: Model,
    stats_by_model: ModelStats,
    model_aggregate: AggregateInfo,
    contributor_lookup: Dict[str, ContributorProviderStats],
    provider_lookup: Dict[str, ContributorProviderStats],
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
        for association in _iter_activity_associations(activity):
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
                                    contributor, ContributorProviderStats()
                                )
                                contributor_info.set_references.add(ref)
                        if provenance.provided_by:
                            for provider in provenance.provided_by:
                                provider_info = provider_lookup.setdefault(
                                    provider, ContributorProviderStats()
                                )
                                provider_info.set_references.add(ref)


def _collect_molecule_terms(
    activity,
    stats_by_model: ModelStats,
    model_aggregate: AggregateInfo,
    contributor_lookup: Dict[str, ContributorProviderStats],
    provider_lookup: Dict[str, ContributorProviderStats],
) -> None:
    """Collect molecule input/output terms from MoleculeAssociation objects.

    Iterates over has_input, has_primary_input, has_output, and
    has_primary_output associations for the given activity, and attributes
    terms to contributors and providers via association-level provenances.

    Populates:
        stats_by_model.list_has_input_term / list_has_output_term
        model_aggregate.list_has_input_term / list_has_output_term
        contributor_lookup[contributor].list_has_input_term / list_has_output_term
        provider_lookup[provider].list_has_input_term / list_has_output_term
    """
    # Collect input molecule associations
    input_associations: list[MoleculeAssociation] = list(activity.has_input or [])
    if activity.has_primary_input:
        input_associations.append(activity.has_primary_input)

    for ma in input_associations:
        if not ma.term:
            continue
        stats_by_model.list_has_input_term.append(ma.term)
        model_aggregate.list_has_input_term.append(ma.term)
        if ma.provenances:
            for provenance in ma.provenances:
                if provenance.contributor:
                    for contributor in provenance.contributor:
                        contributor_info = contributor_lookup.setdefault(
                            contributor, ContributorProviderStats()
                        )
                        contributor_info.list_has_input_term.append(ma.term)
                if provenance.provided_by:
                    for provider in provenance.provided_by:
                        provider_info = provider_lookup.setdefault(
                            provider, ContributorProviderStats()
                        )
                        provider_info.list_has_input_term.append(ma.term)

    # Collect output molecule associations
    output_associations: list[MoleculeAssociation] = list(activity.has_output or [])
    if activity.has_primary_output:
        output_associations.append(activity.has_primary_output)

    for ma in output_associations:
        if not ma.term:
            continue
        stats_by_model.list_has_output_term.append(ma.term)
        model_aggregate.list_has_output_term.append(ma.term)
        if ma.provenances:
            for provenance in ma.provenances:
                if provenance.contributor:
                    for contributor in provenance.contributor:
                        contributor_info = contributor_lookup.setdefault(
                            contributor, ContributorProviderStats()
                        )
                        contributor_info.list_has_output_term.append(ma.term)
                if provenance.provided_by:
                    for provider in provenance.provided_by:
                        provider_info = provider_lookup.setdefault(
                            provider, ContributorProviderStats()
                        )
                        provider_info.list_has_output_term.append(ma.term)


def _collect_terms(
    activity,
    stats_by_model: ModelStats,
    model_aggregate: AggregateInfo,
    contributor_lookup: Dict[str, ContributorProviderStats],
    provider_lookup: Dict[str, ContributorProviderStats],
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

        # Check whether this object carries a GO term
        term = getattr(obj, 'term', None)
        if isinstance(term, str) and term.upper().startswith("GO:"):
            stats_by_model.list_go_terms.append(term)
            model_aggregate.list_go_terms.append(term)
            provenances = getattr(obj, 'provenances', None)
            if provenances:
                for provenance in provenances:
                    if provenance.contributor:
                        for contributor in provenance.contributor:
                            contributor_info = contributor_lookup.setdefault(
                                contributor, ContributorProviderStats()
                            )
                            contributor_info.list_go_terms.append(term)
                    if provenance.provided_by:
                        for provider in provenance.provided_by:
                            provider_info = provider_lookup.setdefault(
                                provider, ContributorProviderStats()
                            )
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


def process_gocam_model_file(
    json_file: Path,
    output_dir: Path | None,
    model_aggregate : AggregateInfo,
    contributor_lookup : Dict[str, ContributorProviderStats],
    provider_lookup: Dict[str, ContributorProviderStats]
) -> ProcessingResult:
    """Process a single GOCAM model JSON file and output statistics information about the model.

    Args:
        json_file (Path): Path to the GO-CAM model JSON file.
        output_dir (Path, optional): Directory to save the GO-CAM model statistics file. If None, no
            file will be written.

    Returns:
        tuple: A tuple indicating the result of processing the GO-CAM file. Possible values are:
            - (ResultType.SUCCESS, None): If the conversion was successful.
            - (ResultType.ERROR, ErrorReason): If there was an error, along with the reason.
    """
    indexer = Indexer()
    # Read GO-CAM JSON file
    try:
        with open(json_file, "r") as f:
            #json_content = f.read()
            #data = json.loads(json_content)
            #model_json = json.dumps(data)
            #gocam_model = Model.model_validate_json(model_json)
            gocam_model = Model.model_validate_json(f.read())
        logger.debug(f"Successfully read GO-CAM model from {json_file}")
    except Exception as e:
        logger.error(f"Error reading file {json_file}", exc_info=e)
        return ResultType.ERROR, None, ErrorReason.READ_ERROR

    # Populate the model with indexing information
    indexer.index_model(gocam_model)

    # Get model statistics information
    calculated_aggregate_values_by_model = gocam_model.query_index

    # Detailed model statistics information
    stats_by_model = ModelStats()
    model_details = ModelDetails()
    model_details.file_name = json_file.name
    model_details.model_id = gocam_model.id
    model_details.model_name = gocam_model.title
    model_aggregate.list_model_details.append(model_details)

    if gocam_model.activities:
        for activity in gocam_model.activities:

            stats_by_model.set_activities.add(activity.id)
            model_aggregate.set_activities.add(activity.id)
            model_details.set_activities.add(activity.id)
            if isinstance(activity.enabled_by, EnabledByGeneProductAssociation):

                stats_by_model.set_activity_unit_gene_product_enablers.add(activity.id)
                model_aggregate.set_activity_unit_gene_product_enablers.add(activity.id)
                if activity.enabled_by.term:
                    stats_by_model.list_enabled_by_gene_product.append(activity.enabled_by.term)
                    stats_by_model.set_enabled_by_gene_product.add(activity.enabled_by.term)
                    stats_by_model.number_of_genes += 1
                    model_aggregate.set_enabled_by_gene_product.add(activity.enabled_by.term)

            elif isinstance(activity.enabled_by, EnabledByProteinComplexAssociation):
                stats_by_model.set_activity_unit_protein_complex_enablers.add(activity.id)
                model_aggregate.set_activity_unit_protein_complex_enablers.add(activity.id)
                if activity.enabled_by.term:
                    stats_by_model.set_protein_complex_in_activity_term.add(activity.enabled_by.term)
                    model_aggregate.set_protein_complex_in_activity_term.add(activity.enabled_by.term)
                if activity.enabled_by.members:
                    for member in activity.enabled_by.members:
                        if member.term:
                            stats_by_model.set_protein_complex_genes.add(member.term)
                            model_aggregate.set_protein_complex_genes.add(member.term)
                            stats_by_model.number_of_genes += 1

            if activity.enabled_by is not None and activity.enabled_by.provenances is not None:
                for provenance in activity.enabled_by.provenances:
                    if provenance.contributor:
                        for contributor in provenance.contributor:
                            contributor_info = contributor_lookup.setdefault(contributor, ContributorProviderStats())
                            _update_entity_evidence_stats(
                                contributor_info, activity, gocam_model.id,
                            )
                    if provenance.provided_by:
                        for provider in provenance.provided_by:
                            provider_info = provider_lookup.setdefault(provider, ContributorProviderStats())
                            _update_entity_evidence_stats(
                                provider_info, activity, gocam_model.id,
                            )

            #Add information about causal relations
            if activity.causal_associations:
                for causal_association in activity.causal_associations:
                    stats_by_model.list_causal_relations.append(activity.id)
                    if causal_association.provenances:
                        for provenance in causal_association.provenances:
                            if provenance.contributor:
                                for contributor in provenance.contributor:
                                    contributor_info = contributor_lookup.setdefault(contributor, ContributorProviderStats())
                                    contributor_info.set_causal_relations.add(activity.id)
                                    contributor_info.number_of_causal_relations += 1
                            if provenance.provided_by:
                                for provider in provenance.provided_by:
                                    provider_info = provider_lookup.setdefault(provider, ContributorProviderStats())
                                    provider_info.set_causal_relations.add(activity.id)
                                    provider_info.number_of_causal_relations += 1
                                    
            # Handle has input and has output
            _collect_molecule_terms(activity, stats_by_model, model_aggregate, contributor_lookup, provider_lookup)

            # Collect GO terms from all associations
            _collect_terms(activity, stats_by_model, model_aggregate, contributor_lookup, provider_lookup)

    # Collect references from all evidence items across all associations
    _collect_references(gocam_model, stats_by_model, model_aggregate, contributor_lookup, provider_lookup)

    #Set fields, counts and sort data for current model
    stats_by_model.number_of_causal_relations = calculated_aggregate_values_by_model.number_of_causal_associations
    stats_by_model.number_of_activitiy_units = calculated_aggregate_values_by_model.number_of_activities
    stats_by_model.number_of_unique_references = len(stats_by_model.set_references)
    stats_by_model.number_of_unique_pmid = _count_pmids(stats_by_model.set_references)

    stats_by_model.number_of_activity_units_enabled_by_gene_product = len(stats_by_model.set_activity_unit_gene_product_enablers)
    stats_by_model.number_of_activity_units_enabled_by_protein_complex = len(stats_by_model.set_activity_unit_protein_complex_enablers)

    stats_by_model.number_of_unique_gene_product_enablers = len(stats_by_model.set_enabled_by_gene_product)
    stats_by_model.number_of_unique_protein_complex_genes = len(stats_by_model.set_protein_complex_genes)
    stats_by_model.number_of_unique_gene_product_and_protein_complex_gene_enablers =  len(stats_by_model.set_enabled_by_gene_product.union(stats_by_model.set_protein_complex_genes))
    stats_by_model.list_enabled_by_gene_product.sort()
    
    stats_by_model.number_of_inputs = len(stats_by_model.list_has_input_term)
    stats_by_model.number_of_chemical_inputs = count_chebis(stats_by_model.list_has_input_term)
    unique_has_input = set(stats_by_model.list_has_input_term)
    stats_by_model.number_of_unique_chemical_inputs = count_chebis(unique_has_input)
    stats_by_model.number_of_other_inputs = stats_by_model.number_of_inputs - stats_by_model.number_of_chemical_inputs
    stats_by_model.number_of_unique_other_inputs = len(unique_has_input) -  stats_by_model.number_of_unique_chemical_inputs
    stats_by_model.number_of_outputs = len(stats_by_model.list_has_output_term)
    stats_by_model.number_of_chemical_outputs = count_chebis(stats_by_model.list_has_output_term)
    unique_has_output = set(stats_by_model.list_has_output_term)
    stats_by_model.number_of_unique_chemical_outputs = count_chebis(unique_has_output)
    stats_by_model.number_of_other_outputs = stats_by_model.number_of_inputs - stats_by_model.number_of_chemical_outputs
    stats_by_model.number_of_unique_other_outputs = len(unique_has_output) -  stats_by_model.number_of_unique_chemical_outputs    
    
    stats_by_model.number_of_go_terms = len(stats_by_model.list_go_terms)
    stats_by_model.number_of_unique_go_terms = len(set(stats_by_model.list_go_terms))


    #Update aggregate model data
    model_aggregate.total_number_of_entities_processed += 1
    model_aggregate.number_of_genes = model_aggregate.number_of_genes + stats_by_model.number_of_genes
    model_aggregate.number_of_causal_relations = model_aggregate.number_of_causal_relations + stats_by_model.number_of_causal_relations

    # If dry run is enabled, skip writing the output file
    if output_dir is None:
        logger.info(
            f"Dry run enabled; skipping write for GO-CAM model {gocam_model.id}"
        )
        return ResultType.SUCCESS, calculated_aggregate_values_by_model, None

    # Write GO-CAM stats to output directory
    calculated_aggregate_values_by_model_subdir = "calculated_aggregate_values_by_model"
    calculated_aggregate_values_by_model_name = "calculated_aggregate_values_by model_" + json_file.name
    calculated_aggregate_values_by_model_name_file = output_dir / calculated_aggregate_values_by_model_subdir/ calculated_aggregate_values_by_model_name
    calculated_aggregate_values_by_model_name_file.parent.mkdir(parents=True, exist_ok=True)

    stats_by_model_subdir = "stats_by_model"
    stats_by_model_name = "stats_by_model_" + json_file.name
    stats_by_model_output_file = output_dir / stats_by_model_subdir/ stats_by_model_name
    stats_by_model_output_file.parent.mkdir(parents=True, exist_ok=True)

    try:
        model_stats_json = calculated_aggregate_values_by_model.model_dump_json(exclude_none=True)
        with open(calculated_aggregate_values_by_model_name_file, "w") as f:
            f.write(model_stats_json)
        logger.info(f"Successfully wrote GO-CAM all_stats to {calculated_aggregate_values_by_model_name_file}")

        stats_by_model_json = stats_by_model.model_dump_json(exclude_none=True)
        with open(stats_by_model_output_file, "w") as f:
            f.write(stats_by_model_json)
        logger.info(f"Successfully wrote GO-CAM dertailed_model_stats to {stats_by_model_output_file}")

    except Exception as e:
        logger.error(f"An exception has occurred: {e}")
        return ResultType.ERROR, None, ErrorReason.WRITE_ERROR

    # If we reach here, the conversion and writing were successful
    return ResultType.SUCCESS, calculated_aggregate_values_by_model, None


def create_filename_from_url(url, replacement='_'):
    """
    Generates a safe filename from a URL.
    """
    # 1. Parse the URL to get the path component
    parsed_url = urlparse(url)
    path = parsed_url.path

    # 2. Extract the base filename from the path
    # Use os.path.basename for handling different path separators, though urlparse usually handles this
    filename = os.path.basename(path)

    # If the URL ends with a '/', filename might be empty, use the hostname as a fallback
    if not filename:
        filename = parsed_url.hostname or 'downloaded_file'

    # 3. Sanitize the filename for cross-platform compatibility
    # Invalid characters are typically /\\:*?"<>| and null character
    safe_filename = re.sub(r'[\\/:*?"<>| \n\t]', replacement, filename)

    # 4. Truncate if necessary (some file systems have length limits, e.g., 255 characters)
    max_length = 250
    if len(safe_filename) > max_length:
        safe_filename = safe_filename[:max_length]

    # Ensure filename is not empty or invalid (e.g., "." or "..")
    if not safe_filename or safe_filename in ('.', '..'):
        safe_filename = 'safe_download'

    return safe_filename

def output_entity_results(
    output_dir: Path | None,
    entity_label: str,
    entity_lookup : Dict[str, ContributorProviderStats],
    entity_sub_dir: str,
    entity_agg_file_name: str,
)-> None:
    """ Output results for entities"""
    entity_agg = AggregateInfo()
    entity_agg.entity = entity_label
    entity_agg.total_number_of_entities_processed = len(entity_lookup)

    for entity, details in entity_lookup.items():
        logger.debug(f"Processing contributor: {entity}")
        #Set numbers
        details.number_of_models = len(details.set_models)
        details.number_of_unique_gene_product_enablers = len(details.set_enabled_by_gene_product)
        details.number_of_unique_member_protein_complex_genes = len(details.set_protein_complex_genes)
        details.number_of_unique_gene_product_and_protein_complex_gene_enablers = len(details.set_enabled_by_gene_product.union(details.set_protein_complex_genes))
        details.number_of_unique_references = len(details.set_references)
        details.number_of_activitiy_units = len(details.set_activities)
        details.number_of_unique_activity_units_enabled_by_gene_product_association = len(details.set_activity_unit_gene_product_enablers)
        details.number_of_unique_activity_units_enabled_by_protein_complex_association = len(details.set_activity_unit_protein_complex_enablers)
        details.number_of_unique_causal_relations = len(details.set_causal_relations)
        details.number_of_inputs = len(details.list_has_input_term)
        details.number_of_chemical_inputs = count_chebis(details.list_has_input_term)
        unique_has_input = set(details.list_has_input_term)
        details.number_of_unique_chemical_inputs = count_chebis(unique_has_input)
        details.number_of_other_inputs = details.number_of_inputs - details.number_of_chemical_inputs
        details.number_of_unique_other_inputs = len(unique_has_input) -  details.number_of_unique_chemical_inputs
        details.number_of_outputs = len(details.list_has_output_term)
        details.number_of_chemical_outputs = count_chebis(details.list_has_output_term)
        unique_has_output = set(details.list_has_output_term)
        details.number_of_unique_chemical_outputs = count_chebis(unique_has_output)
        details.number_of_other_outputs = details.number_of_outputs - details.number_of_chemical_outputs
        details.number_of_unique_other_outputs = len(unique_has_output) -  details.number_of_unique_chemical_outputs    
        
        details.number_of_go_terms = len(details.list_go_terms)
        details.number_of_unique_go_terms = len(set(details.list_go_terms))                

        #Sort lists
        details.list_enabled_by_gene_product.sort()
        details.number_of_unique_pmid = _count_pmids(details.set_references)


        entity_agg.number_of_genes = entity_agg.number_of_genes + len(details.list_enabled_by_gene_product)
        entity_agg.number_of_causal_relations = entity_agg.number_of_causal_relations + details.number_of_causal_relations
        entity_agg.set_activities.update(details.set_activities)
        entity_agg.set_enabled_by_gene_product.update(details.set_enabled_by_gene_product)
        entity_agg.set_references.update(details.set_references)
        entity_agg.set_activity_unit_gene_product_enablers.update(details.set_activity_unit_gene_product_enablers)
        entity_agg.set_activity_unit_protein_complex_enablers.update(details.set_activity_unit_protein_complex_enablers)
        entity_agg.set_protein_complex_in_activity_term.update(details.set_protein_complex_in_activity_term)
        entity_agg.set_protein_complex_genes.update(details.set_protein_complex_genes)
        entity_agg.list_has_input_term.extend(details.list_has_input_term)
        entity_agg.list_has_output_term.extend(details.list_has_output_term)
        entity_agg.list_go_terms.extend(details.list_go_terms)
        
        
        
        
        entity_agg.number_of_inputs = len(entity_agg.list_has_input_term)
        entity_agg.number_of_chemical_inputs = count_chebis(entity_agg.list_has_input_term)
        unique_has_input = set(entity_agg.list_has_input_term)
        entity_agg.number_of_unique_chemical_inputs = count_chebis(unique_has_input)
        entity_agg.number_of_other_inputs = entity_agg.number_of_inputs - entity_agg.number_of_chemical_inputs
        entity_agg.number_of_unique_other_inputs = len(unique_has_input) -  entity_agg.number_of_unique_chemical_inputs           
        entity_agg.number_of_outputs = len(entity_agg.list_has_output_term)
        entity_agg.number_of_chemical_outputs = count_chebis(entity_agg.list_has_output_term)
        unique_has_output = set(entity_agg.list_has_output_term)
        entity_agg.number_of_unique_chemical_outputs = count_chebis(unique_has_output)
        entity_agg.number_of_other_outputs = entity_agg.number_of_outputs - entity_agg.number_of_chemical_outputs
        entity_agg.number_of_unique_other_outputs = len(unique_has_output) -  entity_agg.number_of_unique_chemical_outputs    
        
        entity_agg.number_of_go_terms = len(entity_agg.list_go_terms)
        entity_agg.number_of_unique_go_terms = len(set(entity_agg.list_go_terms))                

        if output_dir is not None:
            stats_by_entity_subdir = entity_sub_dir
            json_file_name = create_filename_from_url(entity +  ".json")
            stats_by_curator_name = entity_sub_dir + "_" + json_file_name
            stats_by_curator_output_file = output_dir / stats_by_entity_subdir/ stats_by_curator_name
            stats_by_curator_output_file.parent.mkdir(parents=True, exist_ok=True)

            try:
                contributor_model_json = details.model_dump_json(exclude_none=True)
                with open(stats_by_curator_output_file, "w") as f:
                    f.write(contributor_model_json)
                logger.info(f"Successfully wrote contributor model to {stats_by_curator_output_file}")
            except Exception as e:
                logger.error(f"Error writing contributor model file {stats_by_curator_output_file}", exc_info=e)
                return ResultType.ERROR, ErrorReason.WRITE_ERROR

    entity_agg.number_of_unique_activity_units = len(entity_agg.set_activities)
    entity_agg.number_of_unique_gene_product_enablers = len(entity_agg.set_enabled_by_gene_product)
    entity_agg.number_of_unique_references = len(entity_agg.set_references)
    entity_agg.number_of_unique_pmid = _count_pmids(entity_agg.set_references)
    entity_agg.number_of_unique_activity_units_enabled_by_protein_complex_association = len(entity_agg.set_activity_unit_protein_complex_enablers)
    entity_agg.number_of_unique_activity_units_enabled_by_gene_product_association = len(entity_agg.set_activity_unit_gene_product_enablers)
    entity_agg.number_of_unique_protein_complex_terms = len(entity_agg.set_protein_complex_in_activity_term)
    entity_agg.number_of_unique_member_protein_complex_genes = len(entity_agg.set_protein_complex_genes)
    entity_agg.number_of_unique_gene_product_and_protein_complex_gene_enablers = len(entity_agg.set_enabled_by_gene_product.union(entity_agg.set_protein_complex_genes))
       

    if entity_agg.total_number_of_entities_processed != 0:
        entity_agg.average_number_of_unique_gene_product_enablers_for_entity = round(entity_agg.number_of_unique_gene_product_enablers /  entity_agg.total_number_of_entities_processed, 2)



    if output_dir is not None:
        aggregate_stats_name = entity_agg_file_name
        aggregate_stats_name_output_file = output_dir / aggregate_stats_name
        aggregate_stats_name_output_file.parent.mkdir(parents=True, exist_ok=True)

        try:
            aggregate_json = entity_agg.model_dump_json(exclude_none=True)
            with open(aggregate_stats_name_output_file, "w") as f:
                f.write(aggregate_json)
            logger.info(f"Successfully wrote entity aggregate stats to {aggregate_stats_name_output_file}")
        except Exception as e:
            logger.error(f"Error writing entity aggregate stats {aggregate_stats_name_output_file}", exc_info=e)
            return ResultType.ERROR, ErrorReason.WRITE_ERROR


def output_summary(
    results: list[tuple[Path, ProcessingResult]],
    output_dir: Path | None,
    model_aggregate: AggregateInfo,
    contributor_lookup : Dict[str, ContributorProviderStats],
    provider_lookup: Dict[str, ContributorProviderStats],
) -> None:
    """Output a summary of the processing results.

    Args:
        results (list): List of tuples containing the JSON file path and processing result.
        contributor_lookup (hashmap): Lookup of curator to annotation details
        provider_lookup (hashmap): Lookup of group to annotation details
    """
    model_aggregate.entity = "Model"
    model_aggregate.number_of_unique_gene_product_enablers = len(model_aggregate.set_enabled_by_gene_product)
    model_aggregate.number_of_unique_references = len(model_aggregate.set_references)
    model_aggregate.number_of_unique_activity_units = len(model_aggregate.set_activities)
    model_aggregate.number_of_unique_activity_units_enabled_by_protein_complex_association = len(model_aggregate.set_activity_unit_protein_complex_enablers)
    model_aggregate.number_of_unique_activity_units_enabled_by_gene_product_association = len(model_aggregate.set_activity_unit_gene_product_enablers)
    model_aggregate.number_of_unique_protein_complex_terms = len(model_aggregate.set_protein_complex_in_activity_term)
    model_aggregate.number_of_unique_member_protein_complex_genes = len(model_aggregate.set_protein_complex_genes)
    model_aggregate.number_of_unique_gene_product_and_protein_complex_gene_enablers = len(model_aggregate.set_enabled_by_gene_product.union(model_aggregate.set_protein_complex_genes))
    model_aggregate.number_of_inputs = len(model_aggregate.list_has_input_term)
    model_aggregate.number_of_chemical_inputs = count_chebis(model_aggregate.list_has_input_term)
    unique_has_input = set(model_aggregate.list_has_input_term)
    model_aggregate.number_of_unique_chemical_inputs = count_chebis(unique_has_input)
    model_aggregate.number_of_other_inputs = model_aggregate.number_of_inputs - model_aggregate.number_of_chemical_inputs
    model_aggregate.number_of_unique_other_inputs = len(unique_has_input) -  model_aggregate.number_of_unique_chemical_inputs
    model_aggregate.number_of_outputs = len(model_aggregate.list_has_output_term)
    model_aggregate.number_of_chemical_outputs = count_chebis(model_aggregate.list_has_output_term)
    unique_has_output = set(model_aggregate.list_has_output_term)
    model_aggregate.number_of_unique_chemical_outputs = count_chebis(unique_has_output)
    model_aggregate.number_of_other_outputs = model_aggregate.number_of_inputs - model_aggregate.number_of_chemical_outputs
    model_aggregate.number_of_unique_other_outputs = len(unique_has_output) -  model_aggregate.number_of_unique_chemical_outputs    
    model_aggregate.number_of_go_terms = len(model_aggregate.list_go_terms)
    model_aggregate.number_of_unique_go_terms = len(set(model_aggregate.list_go_terms))         
    
    
    
    if model_aggregate.total_number_of_entities_processed != 0:
        model_aggregate.average_number_of_unique_gene_product_enablers_for_entity = round(model_aggregate.number_of_unique_gene_product_enablers /  model_aggregate.total_number_of_entities_processed, 2)

    model_aggregate.number_of_unique_pmid = _count_pmids(model_aggregate.set_references)

    if output_dir is not None:
        aggregate_stats_name = "aggregate_model_stats.json"
        aggregate_stats_name_output_file = output_dir / aggregate_stats_name
        aggregate_stats_name_output_file.parent.mkdir(parents=True, exist_ok=True)

        try:
            aggregate_json = model_aggregate.model_dump_json(exclude_none=True)
            with open(aggregate_stats_name_output_file, "w") as f:
                f.write(aggregate_json)
            logger.info(f"Successfully wrote model aggregate stats to {aggregate_stats_name_output_file}")
        except Exception as e:
            logger.error(f"Error writing model aggregate stats {aggregate_stats_name_output_file}", exc_info=e)
            return ResultType.ERROR, ErrorReason.WRITE_ERROR

    #Output information for contributor and provider
    output_entity_results(output_dir=output_dir, entity_label="Curator", entity_lookup=contributor_lookup, entity_sub_dir="stats_by_curator", entity_agg_file_name="aggregate_curator_stats.json")
    output_entity_results(output_dir=output_dir, entity_label="Group", entity_lookup=provider_lookup, entity_sub_dir="stats_by_group", entity_agg_file_name="aggregate_group_stats.json")

    total_count = len(results)
    success_count = sum(1 for _, (result, _, _) in results if result == ResultType.SUCCESS)
    failure_count = total_count - success_count

    tree = Tree(f"[bold]Processed {total_count} models[/bold]")
    if success_count > 0:
        tree.add(f"Successfully output stats for  [b]{success_count}[/b] models", style="green")

    if failure_count > 0:
        tree.add(f"Did not output stats for  [b]{failure_count}[/b] models", style="red")
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
    contributor_lookup : Dict[str, ContributorProviderStats] = {}
    provider_lookup : Dict[str, ContributorProviderStats] = {}

    for json_file in track(
        json_files, description="Processing GO-CAM models and calculating statistics..."
    ):
        logger.debug(f"Processing file: {json_file}")
        result = process_gocam_model_file(json_file, output_dir=output_dir, model_aggregate = model_aggregate, contributor_lookup = contributor_lookup, provider_lookup = provider_lookup)
        results.append((json_file, result))

    # Print result
    output_summary(results, output_dir=output_dir,  model_aggregate = model_aggregate, contributor_lookup = contributor_lookup, provider_lookup = provider_lookup)


if __name__ == "__main__":
    app()