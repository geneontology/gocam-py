"""Utilities for pipeline scripts."""

import io
import json
import logging
from abc import ABC, abstractmethod
from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

import networkx as nx
from rich import print
from rich.logging import RichHandler
from rich.tree import Tree

from gocam.datamodel import Activity, Model, MoleculeAssociation


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


def all_activity_inputs(activity: Activity) -> list[MoleculeAssociation]:
    """Get all inputs for an activity, including primary input if present.

    Args:
        activity: The activity to get inputs for.

    Returns:
        List of all molecule associations that are inputs to the activity.
    """
    inputs = activity.has_input or []
    if activity.has_primary_input:
        inputs.append(activity.has_primary_input)
    return inputs


def all_activity_outputs(activity: Activity) -> list[MoleculeAssociation]:
    """Get all outputs for an activity, including primary output if present.

    Args:
        activity: The activity to get outputs for.

    Returns:
        List of all molecule associations that are outputs of the activity.
    """
    outputs = activity.has_output or []
    if activity.has_primary_output:
        outputs.append(activity.has_primary_output)
    return outputs


def model_to_graph(model: Model) -> nx.DiGraph:
    """Convert a Model to a NetworkX directed graph.

    Nodes are added for activities with an 'enabled_by' property. Edges are added based on:
    - Causal associations between activities.
    - Outputs of one activity serving as inputs to another activity.

    Args:
        model: The model to convert.

    Returns:
        A directed graph representation of the model.
    """
    graph = nx.DiGraph()

    activities_by_input = defaultdict(list)
    activities_by_output = defaultdict(list)
    for activity in model.activities or []:
        for input_ in all_activity_inputs(activity):
            activities_by_input[input_.term].append(activity.id)
        for output in all_activity_outputs(activity):
            activities_by_output[output.term].append(activity.id)

    for activity in model.activities or []:
        if activity.enabled_by is None:
            continue

        graph.add_node(activity.id)

        for causal_assoc in activity.causal_associations or []:
            downstream_activity_id = causal_assoc.downstream_activity
            if downstream_activity_id:
                graph.add_edge(activity.id, downstream_activity_id)

        for output in all_activity_outputs(activity):
            for downstream_activity_id in activities_by_input.get(output.term, []):
                if downstream_activity_id != activity.id:
                    graph.add_edge(activity.id, downstream_activity_id)

    return graph


class FilterReason(str, Enum):
    NO_ACTIVITY_EDGE = "No activity edge"
    USES_COMPLEMENT = "Uses complement"
    NOT_PRODUCTION_MODEL = "Model status is not 'production'"
    NOT_PATHWAY_LIKE = "Model is not pathway-like"


class ErrorReason(str, Enum):
    READ_ERROR = "Read error"
    CONVERSION_ERROR = "Conversion error"
    INDEXING_ERROR = "Indexing error"
    WRITE_ERROR = "Write error"


class PipelineResult(ABC):
    """Base class for pipeline results."""

    @property
    @abstractmethod
    def status(self) -> str:
        """Return the status string for this result."""
        pass

    def get_report_entry(self, model_id: str) -> dict[str, str | list[str]]:
        """Get a report entry for this result.

        Args:
            model_id: The ID of the model associated with this result

        Returns:
            A dictionary representing the report entry.
        """
        return {
            "model_id": model_id,
            "status": self.status,
        }

    def write_to_file(self, file: io.TextIOWrapper, model_id: str) -> None:
        """Write the result as a line of JSON to the given file.

        Args:
            file: The file to write to.
            model_id: The ID of the model associated with this result.
        """
        entry = self.get_report_entry(model_id)
        file.write(json.dumps(entry) + "\n")


@dataclass(frozen=True)
class SuccessResult(PipelineResult):
    """Result for successful processing."""

    data: Any = None
    warnings: list[str] = field(default_factory=list)

    @property
    def status(self) -> str:
        return "success"

    def get_report_entry(self, model_id: str) -> dict[str, str | list[str]]:
        entry = super().get_report_entry(model_id)
        if self.warnings:
            entry["warnings"] = self.warnings
        return entry


@dataclass(frozen=True)
class FilteredResult(PipelineResult):
    """Result for models filtered out during processing."""

    reason: FilterReason

    @property
    def status(self) -> str:
        return "filtered"

    def get_report_entry(self, model_id: str) -> dict[str, str | list[str]]:
        entry = super().get_report_entry(model_id)
        entry["reason"] = self.reason.value
        return entry


@dataclass(frozen=True)
class ErrorResult(PipelineResult):
    """Result for models that failed to process."""

    reason: ErrorReason
    details: str | None = None

    @property
    def status(self) -> str:
        return "error"

    def get_report_entry(self, model_id: str) -> dict[str, str | list[str]]:
        entry = super().get_report_entry(model_id)
        entry["reason"] = self.reason.value
        if self.details:
            entry["details"] = self.details
        return entry


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
