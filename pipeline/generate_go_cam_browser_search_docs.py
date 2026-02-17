#!/usr/bin/env python3
"""
Generate search docs used by the GO-CAM Browser.
"""

import json
import logging
from pathlib import Path
from typing import Annotated

import typer
from _common import (
    ErrorReason,
    ErrorResult,
    PipelineResult,
    ResultSummary,
    SuccessResult,
    get_json_files,
    setup_logger,
)
from glom import Iter, glom
from rich.progress import track

from gocam.datamodel import Model
from gocam.utils import remove_species_code_suffix

app = typer.Typer()

logger = logging.getLogger(__name__)

spec = {
    "id": "id",
    "title": ("title", str.strip),
    "date_modified": "date_modified",
    "status": "status",
    "taxon": "taxon",
    "taxon_label": "query_index.taxon_label",
    "length_of_longest_causal_association_path": "query_index.length_of_longest_causal_association_path",
    "number_of_activities": "query_index.number_of_activities",
    "number_of_strongly_connected_components": "query_index.number_of_strongly_connected_components",
    "enabled_by_gene_labels": (
        "query_index.model_activity_enabled_by_genes",
        Iter("label").filter().map(remove_species_code_suffix).all(),
    ),
    "enabled_by_gene_ids": ("query_index.model_activity_enabled_by_genes", ["id"]),
    "occurs_in_rollup": ("query_index.model_activity_occurs_in_rollup", ["label"]),
    "occurs_in_term_labels": (
        "query_index.model_activity_occurs_in_terms",
        ["label"],
    ),
    "occurs_in_term_ids": ("query_index.model_activity_occurs_in_terms", ["id"]),
    "part_of_rollup": ("query_index.model_activity_part_of_rollup", ["label"]),
    "part_of_term_labels": ("query_index.model_activity_part_of_terms", ["label"]),
    "part_of_term_ids": ("query_index.model_activity_part_of_terms", ["id"]),
    "provided_by_labels": ("query_index.flattened_provided_by", ["label"]),
    "provided_by_ids": ("query_index.flattened_provided_by", ["id"]),
}


def process_gocam_model_file(json_file: Path) -> PipelineResult:
    """Process a GO-CAM model file to generate a search document.

    Args:
        json_file: The path to the GO-CAM model JSON file.

    Returns:
        PipelineResult: A result object indicating success or failure of the processing.
    """
    try:
        with open(json_file, "r") as f:
            indexed_model = Model.model_validate_json(f.read())
    except Exception as e:
        return ErrorResult(
            reason=ErrorReason.READ_ERROR,
            details=str(e),
        )

    try:
        search_doc = glom(indexed_model, spec)
    except Exception as e:
        return ErrorResult(
            reason=ErrorReason.CONVERSION_ERROR,
            details=str(e),
        )
    return SuccessResult(data=search_doc)


@app.command()
def main(
    input_dir: Annotated[
        Path,
        typer.Option(
            exists=True,
            file_okay=False,
            dir_okay=True,
            readable=True,
            help="Directory containing indexed GO-CAM model JSON files.",
        ),
    ],
    output: Annotated[
        typer.FileTextWrite | None,
        typer.Option(
            help="File to write the generated search documents to. Required unless --dry-run is used.",
        ),
    ] = None,
    report_file: Annotated[
        typer.FileTextWrite | None,
        typer.Option(
            help="JSON Lines file to write a detailed report of the processing results.",
        ),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option(
            help="If set, model processing will be performed but no files will be written.",
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
    """Generate search documents for the GO-CAM Browser."""
    setup_logger(verbose)

    # Validate output file if not a dry run
    if not dry_run and output is None:
        raise typer.BadParameter(
            "Output file must be specified unless --dry-run is used."
        )

    # Get list of JSON files in the input directory
    json_files = get_json_files(input_dir, limit=limit)

    result_summary = ResultSummary()
    search_docs = []
    for json_file in track(json_files, description="Generating search docs..."):
        result = process_gocam_model_file(json_file)
        result_summary.add_result(json_file.stem, result)
        if report_file:
            result.write_to_file(report_file, json_file.stem)
        if isinstance(result, SuccessResult):
            search_docs.append(result.data)

    # Print processing summary
    result_summary.print()

    # Write search documents to output file if not a dry run
    if not dry_run and output is not None:
        search_docs.sort(key=lambda m: m.get("date_modified") or "", reverse=True)
        try:
            json.dump(search_docs, output)
        except Exception:
            logger.exception(f"Error writing search documents to file {output.name}")


if __name__ == "__main__":
    app()
