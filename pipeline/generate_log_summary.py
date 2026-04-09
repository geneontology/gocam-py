#!/usr/bin/env python3
"""Generate an Excel summary of pipeline run results based on JSONL log files."""

import itertools
import json
import logging
from collections import defaultdict
from pathlib import Path
from typing import Annotated, Any, Iterable

import typer
from _common import normalize_model_id, setup_logger
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from rich.progress import Progress, track

app = typer.Typer()

logger = logging.getLogger(__name__)


def discover_step_files(logs_dir: Path, extension: str = ".jsonl") -> list[Path]:
    """Find all JSONL files in the logs directory that represent pipeline step reports.

    Args:
        logs_dir: Path to the directory containing log files.
        extension: The file extension to filter log files (default: .jsonl).

    Returns:
        A sorted list of Paths to JSONL files in the logs directory.
    """
    return sorted(
        file
        for file in logs_dir.iterdir()
        if file.is_file() and file.suffix == extension and not file.name.startswith(".")
    )


def iter_log_results(step_file: Path) -> Iterable[tuple[str, dict[str, Any]]]:
    """Iterate over the entries in a JSONL log file, yielding model_id and entry data.

    Args:
        step_file: Path to the JSONL log file.

    Yields:
        Tuples of (normalized_model_id, entry_dict) for each entry in the log file
    """
    with step_file.open() as f:
        for line_number, line in enumerate(f, start=1):
            try:
                entry = json.loads(line)
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"Invalid JSON in {step_file} at line {line_number}: {error}"
                ) from error

            model_id = entry.get("model_id")
            if model_id is None:
                raise ValueError(
                    f"Missing model_id in {step_file} at line {line_number}"
                )
            yield normalize_model_id(model_id), entry


def hyperlink_formula(url: str, display_text: str) -> str:
    """Generate an Excel formula for a hyperlink.

    Args:
        url: The URL the hyperlink should point to.
        display_text: The text to display in the cell.

    Returns:
        A string containing the Excel formula for the hyperlink.
    """
    return f"""=HYPERLINK("{url}", "{display_text}")"""


def format_warning(warning: Any) -> str:
    """Format a warning dictionary into a string for display in Excel.

    Args:
        warning: A dictionary containing warning information, expected to have 'type' and 'message' keys.

    Returns:
        A formatted string representing the warning.
    """
    if isinstance(warning, str):
        return warning
    elif isinstance(warning, dict):
        warning_type = warning.get("type", "")
        message = warning.get("message", "")
        return f"{warning_type}{':' if warning_type else ''}{message}"
    else:
        return str(warning)


@app.command()
def main(
    logs_dir: Annotated[
        Path,
        typer.Option(
            exists=True,
            file_okay=False,
            dir_okay=True,
            readable=True,
            help="Directory containing step JSONL report files.",
        ),
    ],
    output: Annotated[
        Path,
        typer.Option(
            file_okay=True,
            dir_okay=False,
            writable=True,
            help="File to write the generated Excel summary to.",
        ),
    ],
    log_file_extension: Annotated[
        str,
        typer.Option(
            help="File extension used to find log files in logs directory (default: .jsonl).",
        ),
    ] = ".jsonl",
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
            help="Limit the number of models included in the generated summary."
        ),
    ] = 0,
) -> None:
    setup_logger(verbose)

    if not output.name.endswith(".xlsx"):
        raise typer.BadParameter("Output file must have .xlsx extension.")

    step_files = discover_step_files(logs_dir, log_file_extension)
    if not step_files:
        raise typer.BadParameter(
            f"No log files with extension {log_file_extension} found in directory: {logs_dir}"
        )

    step_results_by_model_id: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
    for step_file in step_files:
        for model_id, entry in iter_log_results(step_file):
            step_results_by_model_id[model_id].append(entry)

    # tuple of (column name, column width) for the header row
    columns = [
        ("Model ID", 25),
        ("VPE", 10),
        ("Graph Editor", 15),
        ("Minerva JSON", 15),
        ("Title", 55),
        ("Model Status", 13),
        ("Pipeline Status", 14),
        ("Pipeline Status Details", 35),
        ("Groups", 13),
        ("Longest Path", 14),
        ("Warning Count", 15),
        ("Warnings", 200),
    ]

    wb = Workbook()
    ws = wb.active
    # Write header row
    ws.append([col[0] for col in columns])

    # Set column widths
    for i, (_, width) in enumerate(columns, start=1):
        ws.column_dimensions[get_column_letter(i)].width = width

    # Style values for use later
    fill_green = PatternFill(
        start_color="C6EFCE", end_color="C6EFCE", fill_type="solid"
    )
    alignment_wrapped = Alignment(wrap_text=True, vertical="top")
    font_bold = Font(bold=True)

    model_entries = step_results_by_model_id.items()
    if limit > 0:
        model_entries = itertools.islice(model_entries, limit)

    row = 1
    for model_id, entries in track(
        model_entries, description="Building log summary..."
    ):
        row += 1
        title = None
        model_status = None
        groups = None
        longest_path = None
        pipeline_status = "success"
        pipeline_status_details = None
        warning_count = 0
        warnings = []
        for entry in entries:
            warnings.extend(entry.get("warnings") or [])
            warning_count += len(entry.get("warnings") or [])
            meta = entry.get("meta", {})
            if meta:
                if title is None:
                    title = meta.get("title")
                if model_status is None:
                    model_status = meta.get("status")
                if groups is None:
                    groups = meta.get("groups")
                if longest_path is None:
                    longest_path = meta.get("longest_path")
            if entry.get("status") != "success":
                pipeline_status = entry.get("status", "unknown")
                reason = entry.get("reason")
                if reason:
                    pipeline_status_details = reason
                break
        formatted_groups = ", ".join(groups) if groups else None
        formatted_warnings = (
            "\n".join(format_warning(w) for w in warnings) if warnings else None
        )
        ws.append(
            [
                model_id,
                hyperlink_formula(
                    f"http://noctua.geneontology.org/workbench/noctua-visual-pathway-editor/?model_id=gomodel:{model_id}",
                    "VPE",
                ),
                hyperlink_formula(
                    f"http://noctua.geneontology.org/editor/graph/gomodel:{model_id}",
                    "Graph Editor",
                ),
                hyperlink_formula(
                    f"https://go-public.s3.amazonaws.com/files/go-cam/{model_id}.json",
                    "Minerva JSON",
                ),
                title,
                model_status,
                pipeline_status,
                pipeline_status_details,
                formatted_groups,
                longest_path,
                warning_count,
                formatted_warnings,
            ]
        )
        if pipeline_status == "success":
            for cell in ws[row]:
                cell.fill = fill_green

    with Progress() as progress:
        progress.add_task(description="Writing summary file...", total=None)
        # Add filters to the header row
        ws.auto_filter.ref = f"A1:L{row}"

        # Apply text wrapping and alignment to all cells
        for row in ws.iter_rows():
            for cell in row:
                cell.alignment = alignment_wrapped

        # Make the header row bold
        for cell in ws[1]:
            cell.font = font_bold

        # Save the workbook to the specified output file
        wb.save(output)


if __name__ == "__main__":
    app()
