#!/usr/bin/env python3
"""
Convert Minerva models to GO-CAM models and sort them into True GO-CAMs, Pseudo GO-CAMs, and
excluded models.

This script reads Minerva model JSON files from an input directory and translates each one to
a GO-CAM model. It then checks the translated model, and the warnings produced while
translating it, against the criteria below. Checks run in the order listed, and the first
check a model fails decides what happens to it.

Exclusion checks. A model that fails one of these is not written to any output directory. It
is only recorded in the log file:
  - The model has production status.
  - The model has at least one activity with a causal association to another activity, either
    directly or through input and output associations with the same chemical entity.

Pseudo GO-CAM checks. A model that fails one of these is written to the Pseudo GO-CAM output
directory:
  - The model does not use complement types.
  - The model does not use molecular event individuals.
  - Every activity in the model has a causal association with at least one other activity.
  - The model has evidence for at least one assertion.

A model that passes every check is a True GO-CAM and is written to the output directory.
"""

import dataclasses
import json
import logging
from collections.abc import Iterable
from pathlib import Path
from typing import Annotated

import typer
from _common import (
    ErrorReason,
    ErrorResult,
    FilterCategory,
    FilteredResult,
    FilterReason,
    PipelineLogWriter,
    PipelineResult,
    PipelineStep,
    ResultSummary,
    SuccessResult,
    get_json_files,
    setup_logger,
)
from rich.progress import track

from gocam.datamodel import Model, ModelStateEnum
from gocam.translation import (
    MinervaObject,
    MinervaView,
    MinervaWrapper,
    TranslationWarning,
    WarningType,
)
from gocam.utils import all_evidence, all_provenance, model_to_digraph

app = typer.Typer()

logger = logging.getLogger(__name__)


def extract_minerva_model_metadata(minerva_view: MinervaView) -> dict:
    """Extract model summary metadata from top-level Minerva annotations."""
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
    """Extract complete model summary metadata from a translated GO-CAM model."""
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


def classify_model(
    model: Model, warnings: Iterable[TranslationWarning]
) -> SuccessResult | FilteredResult:
    """Check a translated model against the True GO-CAM criteria.

    The checks and their order are described in the module docstring. Checking stops at the
    first check the model fails.

    Args:
        model: The translated GO-CAM model.
        warnings: The warnings produced while translating the model.

    Returns:
        A SuccessResult if the model passes every check, otherwise a FilteredResult with the
        reason for the first failed check. The result has no metadata or warnings; the caller
        adds them.
    """
    if model.status != ModelStateEnum.production:
        return FilteredResult(reason=FilterReason.NOT_PRODUCTION_MODEL)

    graph = model_to_digraph(model)
    if graph.number_of_edges() == 0:
        return FilteredResult(reason=FilterReason.NO_ACTIVITY_EDGE)

    warning_types = {warning.type for warning in warnings}
    if WarningType.INDIVIDUAL_HAS_COMPLEMENT_TYPE in warning_types:
        return FilteredResult(reason=FilterReason.USES_COMPLEMENT)

    if WarningType.MOLECULAR_EVENT_INDIVIDUAL in warning_types:
        return FilteredResult(reason=FilterReason.USES_MOLECULAR_EVENT)

    if any(degree == 0 for _, degree in graph.degree()):
        return FilteredResult(reason=FilterReason.DISCONNECTED_ACTIVITY)

    if next(all_evidence(model), None) is None:
        return FilteredResult(reason=FilterReason.NO_EVIDENCE)

    return SuccessResult()


def process_minerva_model_file(
    json_file: Path,
    output_dir: Path | None,
    pseudo_gocam_output_dir: Path | None,
) -> PipelineResult:
    """Convert one Minerva model JSON file to a GO-CAM model, check it, and write it out.

    True GO-CAMs are written to output_dir and Pseudo GO-CAMs to pseudo_gocam_output_dir,
    using the input file's name. Excluded models are not written.

    Args:
        json_file: Path to the Minerva model JSON file.
        output_dir: Directory to save True GO-CAM models in. If None, True GO-CAMs are not
            written.
        pseudo_gocam_output_dir: Directory to save Pseudo GO-CAM models in. If None, Pseudo
            GO-CAMs are not written.

    Returns:
        A PipelineResult describing the outcome. It includes the model's metadata whenever
        the metadata could be read from the file.
            - SuccessResult: The model is a True GO-CAM. Includes any translation warnings.
            - FilteredResult: The model failed a check. Includes the reason and any
              translation warnings.
            - ErrorResult: The file could not be read, translated, or written. Includes the
              reason and the error message.
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

    result = dataclasses.replace(
        classify_model(gocam_model, translation_result.warnings),
        warnings=translation_warnings,
        meta=meta,
    )
    if isinstance(result, SuccessResult):
        destination = output_dir
    else:
        logger.info(f"GO-CAM model {gocam_model.id} filtered: {result.reason.value}")
        if result.reason.category == FilterCategory.PSEUDO_GOCAM:
            destination = pseudo_gocam_output_dir
        else:
            destination = None

    # Excluded models are never written, and nothing is written during a dry run
    if destination is None:
        return result

    # Write GO-CAM model to the destination directory
    output_file = destination / json_file.name
    try:
        gocam_model_json = gocam_model.model_dump_json(exclude_none=True)
        with open(output_file, "w") as f:
            f.write(gocam_model_json)
        logger.info(f"Successfully wrote GO-CAM model to {output_file}")
    except Exception as e:
        logger.error(f"Error writing GO-CAM model to file {output_file}", exc_info=e)
        return ErrorResult(reason=ErrorReason.WRITE_ERROR, details=str(e), meta=meta)

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
            help="Directory to save True GO-CAM model files. Required unless --dry-run is used.",
        ),
    ] = None,
    pseudo_gocam_output_dir: Annotated[
        Path | None,
        typer.Option(
            exists=True,
            file_okay=False,
            dir_okay=True,
            writable=True,
            help="Directory to save Pseudo GO-CAM model files. Required unless --dry-run is used.",
        ),
    ] = None,
    log_file: Annotated[
        typer.FileTextWrite | None,
        typer.Option(
            help="JSON Lines log file to write conversion results.",
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
    Convert Minerva models to GO-CAM models. True GO-CAMs are written to the output directory and
    Pseudo GO-CAMs are written to the Pseudo GO-CAM output directory. Excluded models are only
    recorded in the log file.
    """
    setup_logger(verbose)

    # Validate output directories if not a dry run
    if dry_run:
        output_dir = None
        pseudo_gocam_output_dir = None
    else:
        if output_dir is None:
            raise typer.BadParameter(
                "Output directory must be specified unless --dry-run is used."
            )
        if pseudo_gocam_output_dir is None:
            raise typer.BadParameter(
                "Pseudo GO-CAM output directory must be specified unless --dry-run is used."
            )

    log_writer = PipelineLogWriter(log_file, PipelineStep.CONVERT) if log_file else None

    # Get list of JSON files in the input directory
    json_files = get_json_files(input_dir, limit=limit)

    # Process each JSON file
    result_summary = ResultSummary()
    for json_file in track(
        json_files, description="Converting and filtering Minerva models..."
    ):
        logger.debug(f"Processing file: {json_file}")
        result = process_minerva_model_file(
            json_file,
            output_dir=output_dir,
            pseudo_gocam_output_dir=pseudo_gocam_output_dir,
        )
        model_id = json_file.stem
        result_summary.add_result(model_id, result)
        if log_writer:
            log_writer.write_result(result, model_id)

    # Print result
    result_summary.print()


if __name__ == "__main__":
    app()
