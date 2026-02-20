#!/usr/bin/env python3
"""
Filter a collection of models by whether they pass the True GO-CAM model criteria.

A True GO-CAM model is defined as a production GO-CAM model that is pathway-like, meaning it
contains at least three activities connected in sequence. The connection can be either a direct
causal association or an indirect association via shared chemical inputs/outputs.
"""

import logging
import shutil
from pathlib import Path
from typing import Annotated

import networkx as nx
import typer
from _common import (
    ErrorReason,
    ErrorResult,
    FilteredResult,
    FilterReason,
    PipelineResult,
    ResultSummary,
    SuccessResult,
    get_json_files,
    setup_logger,
)
from rich.progress import track

from gocam.datamodel import Model, ModelStateEnum
from gocam.indexing.indexer import model_to_digraph

app = typer.Typer()

logger = logging.getLogger(__name__)


def is_model_pathway_like(model: Model) -> bool:
    """Determine if a model is pathway-like based on its graph representation.

    A model is considered pathway-like if there are at least three activities connected in sequence.

    Args:
        model: The model to evaluate.

    Returns:
        True if the model is pathway-like, False otherwise.
    """

    graph = model_to_digraph(model)
    for node in graph.nodes():
        for other_node in graph.nodes():
            if node != other_node:
                try:
                    shortest_path = nx.shortest_path(
                        graph, source=node, target=other_node
                    )
                except nx.NetworkXNoPath:
                    continue
                if len(shortest_path) >= 3:
                    return True
    return False


def process_gocam_model_file(
    json_file: Path,
    output_dir: Path | None,
    pseudo_gocam_output_dir: Path | None,
) -> PipelineResult:
    """Process a GO-CAM model file to determine if it is a True GO-CAM model.

    Args:
        json_file: Path to the GO-CAM model JSON file.
        output_dir: Directory to save True GO-CAM models.
        pseudo_gocam_output_dir: Directory to save pseudo-GO-CAM models.

    Returns:
        PipelineResult indicating success or filtering reason.
    """
    try:
        with open(json_file, "r") as f:
            model = Model.model_validate_json(f.read())
    except Exception as e:
        logger.error(f"Error reading model from {json_file}: {e}")
        return ErrorResult(reason=ErrorReason.READ_ERROR, details=str(e))

    if model.status != ModelStateEnum.production:
        return FilteredResult(reason=FilterReason.NOT_PRODUCTION_MODEL)

    if is_model_pathway_like(model):
        if output_dir:
            try:
                shutil.copy(json_file, output_dir)
            except Exception as e:
                logger.error(f"Error copying file to {output_dir}: {e}")
                return ErrorResult(reason=ErrorReason.WRITE_ERROR, details=str(e))
        return SuccessResult()
    else:
        if pseudo_gocam_output_dir:
            try:
                shutil.copy(json_file, pseudo_gocam_output_dir)
            except Exception as e:
                logger.error(f"Error copying file to {pseudo_gocam_output_dir}: {e}")
                return ErrorResult(reason=ErrorReason.WRITE_ERROR, details=str(e))
        return FilteredResult(reason=FilterReason.NOT_PATHWAY_LIKE)


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
            help="Directory to save production True GO-CAM model JSON files. Required unless --dry-run is used.",
        ),
    ] = None,
    pseudo_gocam_output_dir: Annotated[
        Path | None,
        typer.Option(
            exists=True,
            file_okay=False,
            dir_okay=True,
            writable=True,
            help="Directory to save production pseudo-GO-CAM model JSON files. Required unless --dry-run is used.",
        ),
    ] = None,
    report_file: Annotated[
        typer.FileTextWrite | None,
        typer.Option(
            help="JSON Lines file to write a detailed report of the identification results.",
        ),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option(
            help="If set, the identification will be performed but no files will be written.",
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
    """
    Filter input models based on production status and pathway-like connectivity. Models that pass
    the True GO-CAM criteria are copied to the specified output directory, while those that do not
    are copied to the pseudo-GO-CAM output directory.
    """
    setup_logger(verbose)

    # Validate output directories if not a dry run
    if not dry_run and output_dir is None:
        raise typer.BadParameter(
            "Output directory must be specified unless --dry-run is used."
        )
    if not dry_run and pseudo_gocam_output_dir is None:
        raise typer.BadParameter(
            "Pseudo-GO-CAM output directory must be specified unless --dry-run is used."
        )

    # Validate report_file name
    if report_file and not report_file.name.endswith(".jsonl"):
        logger.warning(
            "Report file should have a .jsonl extension for JSON Lines format."
        )

    # Get list of JSON files in the input directory
    json_files = get_json_files(input_dir, limit=limit)

    # Process each JSON file
    result_summary = ResultSummary()
    for json_file in track(json_files, description="Filtering True GO-CAM models..."):
        result = process_gocam_model_file(
            json_file, output_dir, pseudo_gocam_output_dir
        )
        result_summary.add_result(json_file.stem, result)
        if report_file:
            result.write_to_file(report_file, json_file.stem)

    # Print result
    result_summary.print()


if __name__ == "__main__":
    app()
