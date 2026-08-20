#!/usr/bin/env python3
"""
Convert Minerva models to GO-CAM models, filtering out models that are definitely not GO-CAM shaped.

This script reads Minerva model JSON files from a specified input directory, converts them to
GO-CAM models, and saves the converted models to an output directory. Models that do not meet very
basic connectivity criteria are filtered out. At least one of the following conditions must be met:

- The model has at least one activity with a direct causal association to another activity.
- The model has at least one activity with an indirect association to another activity via
  input and output associations to the same chemical entity.

The models which pass the connectivity criteria are written to the output directory. They may or may
not be True GO-CAM models. This output set is designed for QC analysis and further downstream
filtering.
"""

import dataclasses
import json
import logging
from pathlib import Path
from typing import Annotated

import typer
from _common import (
    ErrorReason,
    ErrorResult,
    FilteredResult,
    FilterReason,
    PipelineReportWriter,
    PipelineResult,
    PipelineStep,
    ResultSummary,
    SuccessResult,
    get_json_files,
    setup_logger,
)
from rich.progress import track

from gocam.datamodel import Model
from gocam.translation import (
    MinervaObject,
    MinervaView,
    MinervaWrapper,
    WarningType,
)
from gocam.utils import all_provenance, model_to_digraph

app = typer.Typer()

logger = logging.getLogger(__name__)


def extract_minerva_model_metadata(minerva_view: MinervaView) -> dict:
    """Extract report metadata from a Minerva model's top-level annotations."""
    model = minerva_view.raw_json
    annotations = minerva_view.get_annotations(model)
    annotations_multivalued = minerva_view.get_annotations_multivalued(model)

    return {
        "title": annotations.get("title"),
        "status": annotations.get("state"),
        "date_modified": annotations.get("date"),
        "contributors": annotations_multivalued.get("contributor"),
        "provided_by": annotations_multivalued.get("providedBy"),
        "provenance_scope": "model",
    }


def extract_gocam_model_metadata(model: Model) -> dict:
    """Extract complete report metadata from a translated GO-CAM model."""
    contributors = set()
    provided_by = set()
    for provenance in all_provenance(model):
        contributors.update(provenance.contributor or [])
        provided_by.update(provenance.provided_by or [])

    return {
        "title": model.title,
        "status": model.status,
        "date_modified": model.date_modified,
        "contributors": sorted(contributors),
        "provided_by": sorted(provided_by),
        "provenance_scope": "all",
    }


def process_minerva_model_file(
    json_file: Path,
    output_dir: Path | None,
) -> PipelineResult:
    """Process a single Minerva model JSON file and convert it to a GO-CAM model.

    Args:
        json_file: Path to the Minerva model JSON file.
        output_dir: Directory to save the converted GO-CAM model file. If None, no file will be
            written.

    Returns:
        A ProcessingResult indicating the outcome:
            - SuccessResult: If the conversion was successful, with any warnings.
            - FilteredResult: If the model was filtered out, with the reason.
            - ErrorResult: If there was an error, with the reason.
    """
    # Read Minerva model from JSON file
    try:
        with open(json_file, "r") as f:
            minerva_model = json.load(f)
        logger.debug(f"Successfully read Minerva model from {json_file}")
    except Exception as e:
        logger.error(f"Error reading file {json_file}", exc_info=e)
        return ErrorResult(reason=ErrorReason.READ_ERROR, details=str(e))

    meta = None
    try:
        # Extract model-level metadata directly from the Minerva model so it remains
        # available if translation fails.
        minerva_object = MinervaObject.model_validate(minerva_model)
        minerva_view = MinervaView(minerva_object)
        meta = extract_minerva_model_metadata(minerva_view)

        # Convert Minerva model to GO-CAM model
        translation_result = MinervaWrapper.translate(minerva_view)
        gocam_model = translation_result.result
        meta = extract_gocam_model_metadata(gocam_model)
        translation_warnings = [
            dataclasses.asdict(w) for w in translation_result.warnings
        ]
        logger.debug(
            f"Successfully converted Minerva model to GO-CAM model for {json_file}"
        )
    except Exception as e:
        logger.error(
            f"Error converting Minerva model to GO-CAM model for {json_file}",
            exc_info=e,
        )
        return ErrorResult(
            reason=ErrorReason.CONVERSION_ERROR, details=str(e), meta=meta
        )

    # Detect if there is at least one activity edge in the model. If not, skip writing the model.
    graph = model_to_digraph(gocam_model)
    has_activity_edge = graph.number_of_edges() > 0
    if not has_activity_edge:
        logger.info(f"GO-CAM model {gocam_model.id} has no activity edges; skipping.")
        return FilteredResult(reason=FilterReason.NO_ACTIVITY_EDGE, meta=meta)

    # Detect if the Minerva model uses complement constructs. If so, skip writing the model.
    if any(
        w
        for w in translation_result.warnings
        if w.type == WarningType.INDIVIDUAL_HAS_COMPLEMENT_TYPE
    ):
        logger.info(
            f"Minerva model for GO-CAM model {gocam_model.id} uses complement constructs; skipping."
        )
        return FilteredResult(reason=FilterReason.USES_COMPLEMENT, meta=meta)

    # If dry run is enabled, skip writing the output file
    if output_dir is None:
        logger.info(
            f"Dry run enabled; skipping write for GO-CAM model {gocam_model.id}"
        )
        return SuccessResult(warnings=translation_warnings, meta=meta)

    # Write GO-CAM model to output directory
    output_file = output_dir / json_file.name
    try:
        gocam_model_json = gocam_model.model_dump_json(exclude_none=True)
        with open(output_file, "w") as f:
            f.write(gocam_model_json)
        logger.info(f"Successfully wrote GO-CAM model to {output_file}")
    except Exception as e:
        logger.error(f"Error writing GO-CAM model to file {output_file}", exc_info=e)
        return ErrorResult(reason=ErrorReason.WRITE_ERROR, details=str(e))

    # If we reach here, the conversion and writing were successful
    return SuccessResult(warnings=translation_warnings, meta=meta)


@app.command()
def main(
    input_dir: Annotated[
        Path,
        typer.Option(
            exists=True,
            file_okay=False,
            dir_okay=True,
            readable=True,
            help="Directory containing Minerva model JSON files.",
        ),
    ],
    output_dir: Annotated[
        Path | None,
        typer.Option(
            exists=True,
            file_okay=False,
            dir_okay=True,
            writable=True,
            help="Directory to save converted GO-CAM model files. Required unless --dry-run is used.",
        ),
    ] = None,
    report_file: Annotated[
        typer.FileTextWrite | None,
        typer.Option(
            help="JSON Lines file to write a detailed report of the conversion results.",
        ),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option(
            help="If set, the conversion will be performed but no files will be written.",
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
    Convert Minerva models to GO-CAM models, with basic connectivity filtering.
    """
    setup_logger(verbose)

    # Validate output directory if not a dry run
    if not dry_run and output_dir is None:
        raise typer.BadParameter(
            "Output directory must be specified unless --dry-run is used."
        )

    # Validate report_file name
    if report_file and not report_file.name.endswith(".jsonl"):
        logger.warning(
            "Report file should have a .jsonl extension for JSON Lines format."
        )
    report_writer = (
        PipelineReportWriter(report_file, PipelineStep.CONVERT) if report_file else None
    )

    # Get list of JSON files in the input directory
    json_files = get_json_files(input_dir, limit=limit)

    # Process each JSON file
    result_summary = ResultSummary()
    for json_file in track(
        json_files, description="Converting Minerva models to GO-CAM models..."
    ):
        logger.debug(f"Processing file: {json_file}")
        result = process_minerva_model_file(json_file, output_dir=output_dir)
        model_id = json_file.stem
        result_summary.add_result(model_id, result)
        if report_writer:
            report_writer.write_result(result, model_id)

    # Print result
    result_summary.print()


if __name__ == "__main__":
    app()
