#!/usr/bin/env python3
"""
Populate the query_index field of all models in a directory.

This script reads GO-CAM model JSON files from a specified input directory, populates the
query_index field for each model, and saves the updated models to an output directory. The output
directory is not intended to be a release artifact, but the indexed models are used by other
downstream pipeline scripts.
"""

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
from rich.progress import track

from gocam.datamodel import Model
from gocam.indexing.indexer import Indexer

app = typer.Typer()


def process_gocam_model_file(
    json_file: Path, output_dir: Path | None, indexer: Indexer
) -> PipelineResult:
    """Process a GO-CAM model file to populate the query_index field.

    Args:
        json_file: The path to the GO-CAM model JSON file.
        output_dir: The directory to save the updated model file. If None, no file will be written.
        indexer: An instance of the Indexer class used to populate the query_index.

    Returns:
        PipelineResult: A result object indicating success or failure of the processing.
    """
    try:
        with open(json_file, "r") as f:
            model = Model.model_validate_json(f.read())
    except Exception as e:
        return ErrorResult(reason=ErrorReason.READ_ERROR, details=str(e))

    try:
        indexer.index_model(model)
    except Exception as e:
        return ErrorResult(reason=ErrorReason.INDEXING_ERROR, details=str(e))

    if output_dir is not None:
        try:
            with open(output_dir / json_file.name, "w") as f:
                f.write(model.model_dump_json(exclude_none=True))
        except Exception as e:
            return ErrorResult(reason=ErrorReason.WRITE_ERROR, details=str(e))

    return SuccessResult()


@app.command()
def main(
    input_dir: Annotated[
        Path,
        typer.Option(
            exists=True,
            file_okay=False,
            dir_okay=True,
            readable=True,
            help="Directory containing GO-CAM model JSON files.",
        ),
    ],
    output_dir: Annotated[
        Path | None,
        typer.Option(
            exists=True,
            file_okay=False,
            dir_okay=True,
            writable=True,
            help="Directory to save indexed GO-CAM model files. Required unless --dry-run is used.",
        ),
    ] = None,
    report_file: Annotated[
        typer.FileTextWrite | None,
        typer.Option(
            help="JSON Lines file to write a detailed report of the indexing results.",
        ),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option(
            help="If set, the indexing will be performed but no files will be written.",
        ),
    ] = False,
    go_adapter_descriptor: Annotated[
        str,
        typer.Option(
            help="OAK adapter descriptor for GO. See: https://incatools.github.io/ontology-access-kit/packages/selectors.html#ontology-adapter-selectors",
        ),
    ] = "sqlite:obo:go",
    ncbi_taxon_adapter_descriptor: Annotated[
        str,
        typer.Option(
            help="OAK adapter descriptor for the NCBITaxon ontology. See: https://incatools.github.io/ontology-access-kit/packages/selectors.html#ontology-adapter-selectors",
        ),
    ] = "sqlite:obo:ncbitaxon",
    goc_groups_yaml: Annotated[
        Path | None,
        typer.Option(
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            help="YAML file defining GOC groups. If not provided, group information will be fetched from `current.geneontology.org`.",
        ),
    ] = None,
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
) -> None:
    """Populate the query_index field of all models in a directory."""
    setup_logger(verbose)

    # Validate output directories if not a dry run
    if not dry_run and output_dir is None:
        raise typer.BadParameter(
            "Output directory must be specified unless --dry-run is used."
        )

    # Get list of JSON files in the input directory
    json_files = get_json_files(input_dir, limit=limit)

    indexer = Indexer(
        go_adapter_descriptor=go_adapter_descriptor,
        ncbi_taxon_adapter_descriptor=ncbi_taxon_adapter_descriptor,
        goc_groups_yaml_path=goc_groups_yaml,
    )

    # Process each JSON file
    result_summary = ResultSummary()
    for json_file in track(json_files, description="Indexing models..."):
        result = process_gocam_model_file(json_file, output_dir, indexer)
        result_summary.add_result(json_file.stem, result)
        if report_file:
            result.write_to_file(report_file, json_file.stem)

    # Print result summary
    result_summary.print()


if __name__ == "__main__":
    app()
