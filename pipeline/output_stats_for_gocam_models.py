


"""
Output statistics information about GO-CAM models
"""

import json
import logging
from collections import defaultdict
from enum import Enum
from pathlib import Path
from rich import print
from rich.logging import RichHandler
from rich.progress import track
from rich.tree import Tree
#from typing import Annotated, Literal, TypeAlias
from typing import (
    Annotated,
    Any,
    ClassVar,
    Dict,
    List,
    Literal,
    Optional,
    TypeAlias,
    Set,
    Union
)

from linkml_runtime.utils.metamodelcore import (
    bnode,
    empty_dict,
    empty_list
)

from pydantic import (
    Field
)
from gocam.translation import MinervaWrapper
from gocam.datamodel import (
    Activity,
    BiologicalProcessTermObject,
    CellularAnatomicalEntityTermObject,
    CellTypeTermObject,
    ConfiguredBaseModel,
    EvidenceTermObject,
    GeneProductTermObject,
    EnabledByGeneProductAssociation,
    EnabledByProteinComplexAssociation,
    GrossAnatomicalStructureTermObject,
    InformationBiomacromoleculeTermObject,
    Model,
    MolecularFunctionTermObject,
    MoleculeTermObject,
    QueryIndex,
    PhaseTermObject,
    PredicateTermObject,
    ProteinComplexTermObject,
    TaxonTermObject,
    TermObject
)
         


from gocam.indexing.Indexer import Indexer


import typer
import re
from urllib.parse import urlparse

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
    number_of_activities: int = 0
    set_activities: Set[str] = set()

    number_of_unique_genes: int = 0
    set_enabled_by_gene_product: Set[str] = set()
    number_of_genes: int = 0
    list_enabled_by_gene_product: List[str] = [] 

    number_of_unique_references: int = 0
    number_of_unique_pmid: int = 0
    set_references: Set[str] = set()
       
    number_of_causal_relations: int = 0
    list_causal_relations: List[str] = []
    
    set_protein_complex_association_activity: Set[str] = set()
    set_protein_complex_in_activivity_term: Set[str] = set()
    set_protein_complex_genes: Set[str] = set()
        
class ContributorProviderStats(ConfiguredBaseModel):
    number_of_models: int = 0 
    set_models: Set[str] = set()

    number_of_activities: int = 0
    set_activities: Set[str] = set()
    
    number_of_unique_genes: int = 0    
    set_enabled_by_gene_product: Set[str] = set()
    list_enabled_by_gene_product: List[str] = []     
    
    number_of_unique_references: int = 0
    number_of_unique_pmid: int = 0
    set_references: Set[str] = set()
    
    number_of_causal_relations: int = 0
    number_of_unique_causal_relations: int = 0
    set_causal_relations: Set[str] = set()
    
    set_protein_complex_association_activity: Set[str] = set()
    set_protein_complex_in_activivity_term: Set[str] = set()
    set_protein_complex_genes: Set[str] = set()    
                

class AggregateInfo(ConfiguredBaseModel):
    entity: str = ""
    total_number_of_entities_processed: int = 0
    average_number_of_unique_genes_for_entity: float = 0.0
    number_of_genes: int = 0
    number_of_causal_relations: int = 0    
    
    number_of_unique_activities: int = 0
    set_activities: Set[str] = set()
    
    number_of_unique_genes: int = 0
    set_enabled_by_gene_product: Set[str] = set()

    number_of_unique_references: int = 0
    number_of_unique_pmid: int = 0
    set_references: Set[str] = set()
    
    number_of_unique_protein_complex_association_activities: int = 0
    set_protein_complex_association_activity: Set[str] = set()
    
    number_of_unique_protein_complex_terms: int = 0 
    set_protein_complex_in_activivity_term: Set[str] = set()
    
    number_of_unique_member_protein_complex_genes: int = 0    
    set_protein_complex_genes: Set[str] = set()    


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
            json_content = f.read()
            data = json.loads(json_content)
            model_json = json.dumps(data)
            gocam_model = Model.model_validate_json(model_json)
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
    activity_is_enabled_by_gene_product = False
    if gocam_model.activities:
        for activity in gocam_model.activities:
            if isinstance(activity.enabled_by, EnabledByGeneProductAssociation):
                activity_is_enabled_by_gene_product = True
                if activity.enabled_by.term:
                    stats_by_model.list_enabled_by_gene_product.append(activity.enabled_by.term)
                    stats_by_model.set_enabled_by_gene_product.add(activity.enabled_by.term)
                    stats_by_model.number_of_genes += 1
                    model_aggregate.set_enabled_by_gene_product.add(activity.enabled_by.term)
            else:
                activity_is_enabled_by_gene_product = False

            if isinstance(activity.enabled_by, EnabledByProteinComplexAssociation):
                stats_by_model.set_protein_complex_association_activity.add(activity.id)
                model_aggregate.set_protein_complex_association_activity.add(activity.id)
                if activity.enabled_by.term:
                    stats_by_model.set_protein_complex_in_activivity_term.add(activity.enabled_by.term)
                    model_aggregate.set_protein_complex_in_activivity_term.add(activity.enabled_by.term)
                if activity.enabled_by.members:
                    for member in activity.enabled_by.members:
                        stats_by_model.set_protein_complex_genes.add(member)
                        model_aggregate.set_protein_complex_genes.add(member)
            
            stats_by_model.set_activities.add(activity.id)
            model_aggregate.set_activities.add(activity.id)
            
            if activity.enabled_by.evidence:
                for evidence in activity.enabled_by.evidence:
                    if evidence.reference:
                        stats_by_model.set_references.add(evidence.reference)
                        model_aggregate.set_references.add(evidence.reference)
                    if evidence.provenances:
                        for provenance in evidence.provenances:
                            if provenance.contributor:
                                for contributor in provenance.contributor:
                                    if contributor in contributor_lookup:
                                        contributor_info = contributor_lookup[contributor]
                                    else:
                                        contributor_info = ContributorProviderStats()
                                        contributor_lookup[contributor] = contributor_info
                                    
                                    if activity_is_enabled_by_gene_product == True and activity.enabled_by.term:
                                        contributor_info.list_enabled_by_gene_product.append(activity.enabled_by.term)
                                        contributor_info.set_enabled_by_gene_product.add(activity.enabled_by.term)
                                        
                                    if isinstance(activity.enabled_by, EnabledByProteinComplexAssociation):
                                        contributor_info.set_protein_complex_association_activity.add(activity.id)
                                        if activity.enabled_by.term:
                                            contributor_info.set_protein_complex_in_activivity_term.add(activity.enabled_by.term)
                                        if activity.enabled_by.members:
                                            for member in activity.enabled_by.members:
                                                contributor_info.set_protein_complex_genes.add(member)
                        
                                    contributor_info.set_activities.add(activity.id)
                                    contributor_info.set_references.add(evidence.reference)
                                    contributor_info.set_models.add(gocam_model.id)                                            
                                        
                            if provenance.provided_by:
                                for provider in provenance.provided_by:
                                    if provider in provider_lookup:
                                        provider_info = provider_lookup[provider]
                                    else:
                                        provider_info = ContributorProviderStats()
                                        provider_lookup[provider] = provider_info
                                        
                                    if activity_is_enabled_by_gene_product == True and activity.enabled_by.term:
                                        provider_info.list_enabled_by_gene_product.append(activity.enabled_by.term)
                                        provider_info.set_enabled_by_gene_product.add(activity.enabled_by.term)
                                        
                                    if isinstance(activity.enabled_by, EnabledByProteinComplexAssociation):
                                        provider_info.set_protein_complex_association_activity.add(activity.id)
                                        if activity.enabled_by.term:
                                            provider_info.set_protein_complex_in_activivity_term.add(activity.enabled_by.term)
                                        if activity.enabled_by.members:
                                            for member in activity.enabled_by.members:
                                                provider_info.set_protein_complex_genes.add(member)                                        
                                                                                    
                                    provider_info.set_activities.add(activity.id)
                                    provider_info.set_references.add(evidence.reference)
                                    provider_info.set_models.add(gocam_model.id)
            
            #Add information about causal relations
            if activity.causal_associations:
                for causal_association in activity.causal_associations:
                    stats_by_model.list_causal_relations.append(activity.id)
                    if causal_association.provenances:
                        for provenance in causal_association.provenances:        
                            if provenance.contributor:
                                for contributor in provenance.contributor:
                                    if contributor in contributor_lookup:
                                        contributor_info = contributor_lookup[contributor]
                                    else:
                                        contributor_info = ContributorProviderStats()
                                        contributor_lookup[contributor] = contributor_info
                                        
                                    contributor_info.set_causal_relations.add(activity.id)
                                    contributor_info.number_of_causal_relations += 1
                                            
                            if provenance.provided_by:
                                for provider in provenance.provided_by:
                                    if provider in provider_lookup:
                                        provider_info = provider_lookup[provider]
                                    else:
                                        provider_info = ContributorProviderStats()
                                        provider_lookup[provider] = provider_info
                                        
                                    provider_info.set_causal_relations.add(activity.id)
                                    provider_info.number_of_causal_relations += 1  
                            
                            
                        
    #Set fields, counts and sort data for current model
    stats_by_model.number_of_causal_relations = calculated_aggregate_values_by_model.number_of_causal_associations
    stats_by_model.number_of_activities = calculated_aggregate_values_by_model.number_of_activities
    stats_by_model.number_of_unique_references = len(stats_by_model.set_references)
             
    stats_by_model.number_of_unique_genes = len(stats_by_model.set_enabled_by_gene_product)     
    stats_by_model.list_enabled_by_gene_product.sort()     
    
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
    calculated_aggregate_values_by_model_subdir = "calculated_aggregate_values_by model"    
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
    import os
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
    entity_sub_dir,
    entity_agg_file_name,
)-> None:
    """ Output results for entities"""
    entity_agg = AggregateInfo()
    entity_agg.entity = entity_label    
    entity_agg.total_number_of_entities_processed = len(entity_lookup)
    
    for entity, details in entity_lookup.items():
        logger.debug(f"Processing contributor: {entity}")
        #Set numbers
        details.number_of_models = len(details.set_models)
        details.number_of_unique_genes = len(details.set_enabled_by_gene_product)
        details.number_of_unique_references = len(details.set_references)
        details.number_of_activities = len(details.set_activities)
        details.number_of_unique_causal_relations = len(details.set_causal_relations)
        
        #Sort lists
        details.list_enabled_by_gene_product.sort()
        
        for reference in details.set_references:
            if reference.lower().startswith("pmid"):
                details.number_of_unique_pmid += 1
                
                
        entity_agg.number_of_genes = entity_agg.number_of_genes + len(details.list_enabled_by_gene_product)
        entity_agg.number_of_causal_relations = entity_agg.number_of_causal_relations + details.number_of_causal_relations
        entity_agg.set_activities.update(details.set_activities)
        entity_agg.set_enabled_by_gene_product.update(details.set_enabled_by_gene_product)
        entity_agg.set_references.update(details.set_references)
        entity_agg.set_protein_complex_association_activity.update(details.set_protein_complex_association_activity)
        entity_agg.set_protein_complex_in_activivity_term.update(details.set_protein_complex_in_activivity_term)        
        entity_agg.set_protein_complex_genes.update(details.set_protein_complex_genes)    

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
        
    entity_agg.number_of_unique_activities = len(entity_agg.set_activities)
    entity_agg.number_of_unique_genes = len(entity_agg.set_enabled_by_gene_product)
    entity_agg.number_of_unique_references = len(entity_agg.set_references)
    entity_agg.number_of_unique_protein_complex_association_activities = len(entity_agg.set_protein_complex_association_activity)
    entity_agg.number_of_unique_protein_complex_terms = len(entity_agg.set_protein_complex_in_activivity_term)    
    entity_agg.number_of_unique_member_protein_complex_genes = len(entity_agg.set_protein_complex_genes)        
    
    if entity_agg.total_number_of_entities_processed != 0:
        entity_agg.average_number_of_unique_genes_for_entity = round(entity_agg.number_of_unique_genes /  entity_agg.total_number_of_entities_processed, 2)
    
    
    if entity_agg.set_references:
        for reference in entity_agg.set_references:
            if reference.lower().startswith("pmid"):
                entity_agg.number_of_unique_pmid += 1
    
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
    max_ids: int = 5,    
) -> None:
    """Output a summary of the processing results.

    Args:
        results (list): List of tuples containing the JSON file path and processing result.
        max_ids (int): Maximum number of model IDs to display per reason.
        contributor_lookup (hashmap): Lookup of curator to annotation details
        provider_lookup (hashmap): Lookup of group to annotation details
    """
    model_aggregate.entity = "Model"
    model_aggregate.number_of_unique_genes = len(model_aggregate.set_enabled_by_gene_product)
    model_aggregate.number_of_unique_references = len(model_aggregate.set_references)
    model_aggregate.number_of_unique_activities = len(model_aggregate.set_activities)
    model_aggregate.number_of_unique_protein_complex_association_activities = len(model_aggregate.set_protein_complex_association_activity)
    model_aggregate.number_of_unique_protein_complex_terms = len(model_aggregate.set_protein_complex_in_activivity_term)
    model_aggregate.number_of_unique_member_protein_complex_genes = len(model_aggregate.set_protein_complex_genes)        
    if model_aggregate.total_number_of_entities_processed != 0:
        model_aggregate.average_number_of_unique_genes_for_entity = round(model_aggregate.number_of_unique_genes /  model_aggregate.total_number_of_entities_processed, 2)
    

    for reference in model_aggregate.set_references:
        if reference.lower().startswith("pmid"):
            model_aggregate.number_of_unique_pmid += 1
                
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
    error_by_reason: defaultdict[ErrorReason, list[str]] = defaultdict(list)
    for json_file, (result, model_stats, reason) in results:
        model_id = json_file.stem


    error_count = sum(len(v) for v in error_by_reason.values())

    tree = Tree(f"[bold]Processed {total_count} models[/bold]")
    if success_count > 0:
        tree.add(f"Successfully output stats for  [b]{success_count}[/b] models", style="green")



    if error_count > 0:
        error_branch = tree.add(
            f"Failed to output statistics for  [b]{error_count}[/b] models for the following reasons:",
            style="red",
        )
        for reason, ids in error_by_reason.items():
            reason_branch = error_branch.add(
                f"{reason.value}: [b]{len(ids)}[/b] models"
            )
            for model_id in ids[:max_ids]:
                reason_branch.add(model_id)
            if len(ids) > max_ids:
                reason_branch.add(f"... and {len(ids) - max_ids} more.")

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