#!/usr/bin/env python3
"""
Filter a collection of models by whether they pass the True GO-CAM model criteria.

A True GO-CAM model is defined as a model where all of the following criteria are met:
  - The model is in production status.
  - The model has at least two activities that are connected by a causal association, either
    directly or indirectly via shared chemical entities.
  - The model has no activities that are disconnected from all other activities in the model.
"""

import logging
import shutil
from pathlib import Path
from typing import Annotated

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
from gocam.utils import model_to_digraph

app = typer.Typer()

logger = logging.getLogger(__name__)


def is_true_gocam(model: Model) -> PipelineResult:
    """Determine if a GO-CAM model meets the criteria for being a True GO-CAM model.

    Args:
        model: The model to evaluate.

    Returns:
        PipelineResult: A SuccessResult if the model is a True GO-CAM model, or a FilteredResult
        with the reason for filtering if it is not.
    """

    # Check if model is in production status
    if model.status != ModelStateEnum.production:
        return FilteredResult(reason=FilterReason.NOT_PRODUCTION_MODEL)

    graph = model_to_digraph(model)

    # Check for at least one edge (causal association) between activities
    if graph.number_of_edges() < 1:
        return FilteredResult(reason=FilterReason.NO_ACTIVITY_EDGE)

    # Check for nodes with degree 0 (disconnected activities)
    for _, degree in graph.degree():
        if degree == 0:
            return FilteredResult(reason=FilterReason.DISCONNECTED_ACTIVITY)

    return SuccessResult()


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

    result = is_true_gocam(model)
    if isinstance(result, SuccessResult):
        if output_dir:
            try:
                shutil.copy(json_file, output_dir)
            except Exception as e:
                logger.error(f"Error copying file to {output_dir}: {e}")
                return ErrorResult(reason=ErrorReason.WRITE_ERROR, details=str(e))
    elif isinstance(result, FilteredResult):
        if (
            pseudo_gocam_output_dir
            and result.reason != FilterReason.NOT_PRODUCTION_MODEL
        ):
            try:
                shutil.copy(json_file, pseudo_gocam_output_dir)
            except Exception as e:
                logger.error(f"Error copying file to {pseudo_gocam_output_dir}: {e}")
                return ErrorResult(reason=ErrorReason.WRITE_ERROR, details=str(e))
    return result


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
    Filter input models based on True GO-CAM criteria. Models that pass the True GO-CAM criteria are
    copied to the specified output directory, while those that do not are copied to the
    pseudo-GO-CAM output directory. Models that do not have production status are not copied to
    either output directory.
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
