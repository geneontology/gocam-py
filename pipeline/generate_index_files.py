#!/usr/bin/env python3
"""
Generate various index files (JSON) from a directory of GO-CAM models.
"""

import json
import logging
from abc import ABC, abstractmethod
from collections import defaultdict
from pathlib import Path
from typing import Annotated

import typer
from _common import get_json_files, setup_logger
from rich.progress import track

from gocam.datamodel import Model
from gocam.indexing.Indexer import Indexer

app = typer.Typer()

logger = logging.getLogger(__name__)


class EmptyIndexError(Exception):
    """Custom exception raised when an index is empty and will not be written to a file."""


class IndexReport(ABC):
    """
    Base class for generating an index from GO-CAM models.
    """

    def __init__(self, file_name: str):
        self.file_name = file_name
        self.index: dict[str, set] = defaultdict(set)

    @abstractmethod
    def process_model(self, model: Model) -> None:
        """Process a GO-CAM model and update the index accordingly.

        Args:
            model: The GO-CAM model to process.
        """
        pass

    def add(self, key: str, model: Model) -> None:
        """Add a model to the index under the specified key.

        Args:
            key: The key to index the model under.
            model: The GO-CAM model to add to the index.
        """
        model_id = model.id.removeprefix("gomodel:")
        self.index[key].add(model_id)

    def write(self, output_directory: Path) -> None:
        """Write the index to a JSON file.

        Args:
            output_directory: The directory where the index file should be saved.

        Raises:
            EmptyIndexError: If the index is empty
        """
        if not self.index:
            raise EmptyIndexError(
                f"Index for {self.file_name} is empty. No file will be written."
            )

        with open(output_directory / self.file_name, "w") as f:
            json.dump({k: sorted(v) for k, v in self.index.items()}, f, indent=2)


class ContributorIndexReport(IndexReport):
    """Index report for the 'contributor' field in GO-CAM models."""

    def process_model(self, model: Model) -> None:
        if model.query_index is None:
            return

        for contributor in model.query_index.flattened_contributors or []:
            self.add(contributor, model)


class EntityIndexReport(IndexReport):
    """Index report for entities in GO-CAM models, including ancestors of GO terms."""

    def process_model(self, model: Model) -> None:
        if model.query_index is None:
            return

        for entity in model.objects or []:
            self.add(entity.id, model)

        for term in model.query_index.model_activity_molecular_function_closure or []:
            self.add(term.id, model)
        for term in model.query_index.model_activity_part_of_closure or []:
            self.add(term.id, model)
        for term in model.query_index.model_activity_occurs_in_closure or []:
            self.add(term.id, model)


class EvidenceIndexReport(IndexReport):
    """Index report for the 'evidence' field in GO-CAM models."""

    def process_model(self, model: Model) -> None:
        if model.query_index is None:
            return

        for evidence in model.query_index.flattened_evidence_terms or []:
            self.add(evidence.id, model)


class ProvidedByIndexReport(IndexReport):
    """Index report for the 'provided_by' field in GO-CAM models."""

    def process_model(self, model: Model) -> None:
        if model.query_index is None:
            return

        for provided_by in model.query_index.flattened_provided_by or []:
            self.add(provided_by.id, model)


class SourceIndexReport(IndexReport):
    """Index report for the 'references' field in GO-CAM models."""

    def process_model(self, model: Model) -> None:
        if model.query_index is None:
            return

        for ref in model.query_index.flattened_references or []:
            self.add(ref.id, model)


class TaxonIndexReport(IndexReport):
    """Index report for the 'taxon' and 'additional_taxa' fields in GO-CAM models."""

    def process_model(self, model: Model) -> None:
        if model.taxon:
            self.add(model.taxon, model)
        for taxon in model.additional_taxa or []:
            self.add(taxon, model)


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
            help="Directory to save generated index files. Required unless --dry-run is used.",
        ),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option(
            help="If set, model processing will be performed but no files will be written.",
        ),
    ] = False,
    go_adapater_descriptor: Annotated[
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
):
    """
    Generate various indexes (contributor, entity, provided_by, source, taxon) for a directory
    of GO-CAM models.
    """
    setup_logger(verbose)

    # Validate output directories if not a dry run
    if not dry_run and output_dir is None:
        raise typer.BadParameter(
            "Output directory must be specified unless --dry-run is used."
        )

    # Get list of JSON files in the input directory
    json_files = get_json_files(input_dir, limit=limit)

    indexer = Indexer(
        go_adapter_descriptor=go_adapater_descriptor,
        ncbi_taxon_adapter_descriptor=ncbi_taxon_adapter_descriptor,
        goc_groups_yaml_path=goc_groups_yaml,
    )

    reports = [
        ContributorIndexReport("contributor_index.json"),
        EntityIndexReport("entity_index.json"),
        EvidenceIndexReport("evidence_index.json"),
        ProvidedByIndexReport("provided_by_index.json"),
        SourceIndexReport("source_index.json"),
        TaxonIndexReport("taxon_index.json"),
    ]

    for json_file in track(json_files, description="Generating index files..."):
        try:
            with open(json_file, "r") as f:
                model = Model.model_validate_json(f.read())
        except Exception:
            logger.exception(f"Error reading model from {json_file}")
            continue

        indexer.index_model(model)
        if model.query_index is None:
            logger.error(f"Model {json_file} has no query index after indexing.")
            continue

        for report in reports:
            try:
                report.process_model(model)
            except Exception:
                logger.exception(
                    f"Error processing model {json_file} for report {report.file_name}"
                )

    if output_dir and not dry_run:
        for report in reports:
            try:
                report.write(output_dir)
            except Exception:
                logger.exception(f"Error writing index file {report.file_name}")


if __name__ == "__main__":
    app()
