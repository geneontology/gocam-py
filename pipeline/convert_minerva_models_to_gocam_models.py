#!/usr/bin/env python3
"""
Minerva Model to GO-CAM Model Converter
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
from typing import Annotated, Literal, TypeAlias

from gocam.translation import MinervaWrapper

import typer

app = typer.Typer()

logger = logging.getLogger(__name__)


class ResultType(str, Enum):
    SUCCESS = "Success"
    FILTERED = "Filtered"
    ERROR = "Error"


class FilterReason(str, Enum):
    NO_ACTIVITY_EDGE = "No activity edge"
    USES_COMPLEMENT = "Uses complement"


class ErrorReason(str, Enum):
    READ_ERROR = "Read error"
    CONVERSION_ERROR = "Conversion error"
    WRITE_ERROR = "Write error"


ProcessingResult: TypeAlias = (
    tuple[Literal[ResultType.SUCCESS], None]
    | tuple[Literal[ResultType.FILTERED], FilterReason]
    | tuple[Literal[ResultType.ERROR], ErrorReason]
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


def minerva_model_uses_complement(minerva_model: dict) -> bool:
    """Check if the Minerva model uses complement constructs.

    Args:
        minerva_model (dict): The Minerva model as a dictionary.

    Returns:
        bool: True if the model uses complement constructs, False otherwise.
    """
    for individual in minerva_model.get("individuals", []):
        for type_ in individual.get("type", []):
            if type_.get("type") == "complement":
                return True
    return False


def process_minerva_model_file(
    json_file: Path,
    output_dir: Path | None,
) -> ProcessingResult:
    """Process a single Minerva model JSON file and convert it to a GO-CAM model.

    Args:
        json_file (Path): Path to the Minerva model JSON file.
        output_dir (Path, optional): Directory to save the converted GO-CAM model file. If None, no
            file will be written.

    Returns:
        tuple: A tuple indicating the result of the processing. Possible values are:
            - (ResultType.SUCCESS, None): If the conversion was successful.
            - (ResultType.FILTERED, FilterReason): If the model was filtered out, along with the reason.
            - (ResultType.ERROR, ErrorReason): If there was an error, along with the reason.
    """
    # Read Minerva model from JSON file
    try:
        with open(json_file, "r") as f:
            minerva_model = json.load(f)
        logger.debug(f"Successfully read Minerva model from {json_file}")
    except Exception as e:
        logger.error(f"Error reading file {json_file}", exc_info=e)
        return ResultType.ERROR, ErrorReason.READ_ERROR

    # Convert Minerva model to GO-CAM model
    try:
        gocam_model = MinervaWrapper.minerva_object_to_model(minerva_model)
        logger.debug(
            f"Successfully converted Minerva model to GO-CAM model for {json_file}"
        )
    except Exception as e:
        logger.error(
            f"Error converting Minerva model to GO-CAM model for {json_file}", exc_info=e
        )
        return ResultType.ERROR, ErrorReason.CONVERSION_ERROR

    # Detect if there is at least one activity edge in the model. If not, skip writing the model.
    has_activity_edge = False
    for activity in gocam_model.activities or []:
        if activity.causal_associations:
            has_activity_edge = True
            break
    if not has_activity_edge:
        logger.info(f"GO-CAM model {gocam_model.id} has no activity edges; skipping.")
        return ResultType.FILTERED, FilterReason.NO_ACTIVITY_EDGE

    # Detect if the Minerva model uses complement constructs. If so, skip writing the model.
    if minerva_model_uses_complement(minerva_model):
        logger.info(
            f"Minerva model for GO-CAM model {gocam_model.id} uses complement constructs; skipping."
        )
        return ResultType.FILTERED, FilterReason.USES_COMPLEMENT

    # If dry run is enabled, skip writing the output file
    if output_dir is None:
        logger.info(
            f"Dry run enabled; skipping write for GO-CAM model {gocam_model.id}"
        )
        return ResultType.SUCCESS, None

    # Write GO-CAM model to output directory
    output_file = output_dir / json_file.name
    try:
        gocam_model_json = gocam_model.model_dump_json(exclude_none=True)
        with open(output_file, "w") as f:
            f.write(gocam_model_json)
        logger.info(f"Successfully wrote GO-CAM model to {output_file}")
    except Exception as e:
        logger.error(f"Error writing GO-CAM model to file {output_file}", exc_info=e)
        return ResultType.ERROR, ErrorReason.WRITE_ERROR

    # If we reach here, the conversion and writing were successful
    return ResultType.SUCCESS, None


def print_summary(
    results: list[tuple[Path, ProcessingResult]], max_ids: int = 5
) -> None:
    """Print a summary of the processing results.

    Args:
        results (list): List of tuples containing the JSON file path and processing result.
        max_ids (int): Maximum number of model IDs to display per reason.
    """
    total_count = len(results)
    success_count = sum(1 for _, (result, _) in results if result == ResultType.SUCCESS)
    filtered_by_reason: defaultdict[FilterReason, list[str]] = defaultdict(list)
    error_by_reason: defaultdict[ErrorReason, list[str]] = defaultdict(list)
    for json_file, (result, reason) in results:
        model_id = json_file.stem
        if result == ResultType.FILTERED:
            filtered_by_reason[reason].append(model_id)
        elif result == ResultType.ERROR:
            error_by_reason[reason].append(model_id)
    filtered_count = sum(len(v) for v in filtered_by_reason.values())
    error_count = sum(len(v) for v in error_by_reason.values())

    tree = Tree(f"[bold]Processed {total_count} models[/bold]")
    if success_count > 0:
        tree.add(f"Successfully converted [b]{success_count}[/b] models", style="green")

    if filtered_count > 0:
        filtered_branch = tree.add(
            f"Filtered out [b]{filtered_count}[/b] models for the following reasons:",
            style="yellow",
        )
        for reason, ids in filtered_by_reason.items():
            reason_branch = filtered_branch.add(
                f"{reason.value}: [b]{len(ids)}[/b] models"
            )
            for model_id in ids[:max_ids]:
                reason_branch.add(model_id)
            if len(ids) > max_ids:
                reason_branch.add(
                    f"... and {len(ids) - max_ids} more. Run with `--verbose` option and inspect logs to see all."
                )

    if error_count > 0:
        error_branch = tree.add(
            f"Failed to convert [b]{error_count}[/b] models for the following reasons:",
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
    Convert Minerva models to GO-CAM models
    """
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
    for json_file in track(
        json_files, description="Converting Minerva models to GO-CAM models..."
    ):
        logger.debug(f"Processing file: {json_file}")
        result = process_minerva_model_file(json_file, output_dir=output_dir)
        results.append((json_file, result))

    # Print result
    print_summary(results)


if __name__ == "__main__":
    app()
