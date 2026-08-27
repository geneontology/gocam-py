"""Utilities for pipeline scripts."""

import json
import logging
from abc import ABC, abstractmethod
from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, TextIO

from rich import print
from rich.logging import RichHandler
from rich.tree import Tree

logger = logging.getLogger(__name__)


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


def get_json_files(
    input_dir: Path, *, limit: int | None = None, raise_on_empty: bool = True
) -> list[Path]:
    """Get a sorted list of JSON files in the input directory, excluding hidden files.

    Args:
        input_dir: Directory to search for JSON files.
        limit: Optional limit on the number of files to return. If None or 0, no limit is applied.
        raise_on_empty: If True, raise a ValueError if no JSON files are found.

    Returns:
        List of JSON file paths.
    """
    files = sorted(
        file
        for file in input_dir.iterdir()
        if file.is_file() and file.suffix == ".json" and not file.name.startswith(".")
    )
    if not files and raise_on_empty:
        raise ValueError(f"No JSON files found in directory: {input_dir}")
    if limit is not None and limit > 0:
        files = files[:limit]
    return files


def normalize_model_id(model_id: str) -> str:
    """Normalize a model ID by removing the "gomodel:" prefix if present.

    Args:
        model_id: The model ID to normalize.

    Returns:
        The normalized model ID without the "gomodel:" prefix.
    """
    return model_id.removeprefix("gomodel:")


class FilterReason(str, Enum):
    DISCONNECTED_ACTIVITY = "Model has at least one activity with no causal relations"
    NOT_PRODUCTION_MODEL = "Model status is not 'production'"
    NO_ACTIVITY_EDGE = "Model has 0 connected activities"
    NO_EVIDENCE = "Model has no evidence for any assertions"
    USES_COMPLEMENT = "Model uses complement"


class ErrorReason(str, Enum):
    READ_ERROR = "Read error"
    CONVERSION_ERROR = "Conversion error"
    INDEXING_ERROR = "Indexing error"
    WRITE_ERROR = "Write error"


class PipelineStep(str, Enum):
    """Identifiers for steps that produce pipeline logs."""

    CONVERT = "convert"
    FILTER = "filter"
    QUERY_INDEX = "query-index"
    INDEX_FILES = "index-files"
    BROWSER_SEARCH = "browser-search"


PIPELINE_STEP_ORDER = tuple(PipelineStep)
PIPELINE_LOG_RECORD_TYPE = "pipeline_step_log"
PIPELINE_LOG_FORMAT_VERSION = 1


@dataclass(frozen=True, kw_only=True)
class PipelineResult(ABC):
    """Base class for pipeline results."""

    meta: dict[str, Any] | None = None

    @property
    @abstractmethod
    def status(self) -> str:
        """Return the status string for this result."""
        pass

    def get_log_entry(self, model_id: str) -> dict[str, Any]:
        """Get a log entry for this result.

        Args:
            model_id: The ID of the model associated with this result

        Returns:
            A dictionary representing the log entry.
        """
        entry: dict[str, Any] = {
            "model_id": model_id,
            "status": self.status,
        }
        if self.meta:
            entry["meta"] = self.meta
        return entry


@dataclass(frozen=True, kw_only=True)
class SuccessResult(PipelineResult):
    """Result for successful processing."""

    data: Any = None
    warnings: list[Any] = field(default_factory=list)

    @property
    def status(self) -> str:
        return "success"

    def get_log_entry(self, model_id: str) -> dict[str, str | list[str]]:
        entry = super().get_log_entry(model_id)
        if self.warnings:
            entry["warnings"] = self.warnings
        return entry


@dataclass(frozen=True, kw_only=True)
class FilteredResult(PipelineResult):
    """Result for models filtered out during processing."""

    reason: FilterReason

    @property
    def status(self) -> str:
        return "filtered"

    def get_log_entry(self, model_id: str) -> dict[str, str | list[str]]:
        entry = super().get_log_entry(model_id)
        entry["reason"] = self.reason.value
        return entry


@dataclass(frozen=True, kw_only=True)
class ErrorResult(PipelineResult):
    """Result for models that failed to process."""

    reason: ErrorReason
    details: str | None = None

    @property
    def status(self) -> str:
        return "error"

    def get_log_entry(self, model_id: str) -> dict[str, str | list[str]]:
        entry = super().get_log_entry(model_id)
        entry["reason"] = self.reason.value
        if self.details:
            entry["details"] = self.details
        return entry


class PipelineLogWriter:
    """Writer for pipeline step logs in JSON Lines format with identifying header."""

    def __init__(self, file: TextIO, step: PipelineStep) -> None:
        if not file.name.endswith(".jsonl"):
            logger.warning(
                "Log file should have a .jsonl extension for JSON Lines format."
            )
        self.file = file
        self.file.write(
            json.dumps(
                {
                    "record_type": PIPELINE_LOG_RECORD_TYPE,
                    "format_version": PIPELINE_LOG_FORMAT_VERSION,
                    "step": step.value,
                }
            )
            + "\n"
        )

    def write_result(self, result: PipelineResult, model_id: str) -> None:
        """Append one model result to the log."""
        self.file.write(json.dumps(result.get_log_entry(model_id)) + "\n")


class ResultSummary:
    """Helper class to track and summarize pipeline results."""

    def __init__(self, max_ids: int = 5):
        """Initialize the result summary.

        Args:
            max_ids: Maximum number of model IDs to display per reason in the summary.
        """
        self._max_ids = max_ids
        self._total_count = 0
        self._success_count = 0
        self._success_with_warnings_count = 0
        self._filtered_count = 0
        self._error_count = 0
        self._filtered_count_by_reason: defaultdict[str, int] = defaultdict(int)
        self._filtered_ids_by_reason: defaultdict[str, list[str]] = defaultdict(list)
        self._error_count_by_reason: defaultdict[str, int] = defaultdict(int)
        self._error_ids_by_reason: defaultdict[str, list[str]] = defaultdict(list)

    def add_result(self, model_id: str, result: PipelineResult) -> None:
        """Add a processing result to the summary.

        Args:
            model_id: The ID of the model that was processed.
            result: The processing result.
        """
        self._total_count += 1
        match result:
            case SuccessResult(warnings=warnings):
                self._success_count += 1
                if warnings:
                    self._success_with_warnings_count += 1
            case FilteredResult(reason=reason):
                self._filtered_count += 1
                self._filtered_count_by_reason[reason.value] += 1
                if len(self._filtered_ids_by_reason[reason.value]) < self._max_ids:
                    self._filtered_ids_by_reason[reason.value].append(model_id)
            case ErrorResult(reason=reason):
                self._error_count += 1
                self._error_count_by_reason[reason.value] += 1
                if len(self._error_ids_by_reason[reason.value]) < self._max_ids:
                    self._error_ids_by_reason[reason.value].append(model_id)

    def print(self):
        """Print the result summary as a tree to the console."""

        def handle_reasons(
            parent_node: Tree, counts: dict[str, int], ids: dict[str, list[str]]
        ) -> None:
            for reason, total in counts.items():
                reason_branch = parent_node.add(f"{reason}: [b]{total}[/b] models")
                id_list = ids.get(reason, [])
                for model_id in id_list:
                    reason_branch.add(model_id)
                if total > self._max_ids:
                    reason_branch.add(f"... and {total - self._max_ids} more.")

        tree = Tree(f"[bold]Processed {self._total_count} models[/bold]")
        if self._success_count > 0:
            success_branch = tree.add(
                f"Successfully processed [b]{self._success_count}[/b] models",
                style="green",
            )
            if self._success_with_warnings_count > 0:
                success_branch.add(
                    f"With warnings: [b]{self._success_with_warnings_count}[/b] models"
                )

        if self._filtered_count > 0:
            filtered_branch = tree.add(
                f"Filtered out [b]{self._filtered_count}[/b] models for the following reasons:",
                style="yellow",
            )
            handle_reasons(
                filtered_branch,
                self._filtered_count_by_reason,
                self._filtered_ids_by_reason,
            )

        if self._error_count > 0:
            error_branch = tree.add(
                f"Failed to process [b]{self._error_count}[/b] models for the following reasons:",
                style="red",
            )
            handle_reasons(
                error_branch, self._error_count_by_reason, self._error_ids_by_reason
            )

        print(tree)
