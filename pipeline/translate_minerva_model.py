#!/usr/bin/env python3
"""
Translate a single Minerva model JSON file to a GO-CAM model, for debugging translation issues.

The translated GO-CAM model is written to stdout as JSON. Translation warnings are written to
stderr, one per line, as tab-separated values: warning type, entity ID, and message.

This script is a development aid only. It intentionally does not represent a pipeline step: it
does no filtering, writes no pipeline log, and is not run by `run_pipeline.sh` or the GO Release
pipeline. Pipeline steps should not call or depend on it.
"""

import json
import sys
from pathlib import Path
from typing import Annotated

import typer

from gocam.translation import MinervaWrapper, TranslationWarning

app = typer.Typer()


def format_warning(warning: TranslationWarning) -> str:
    """Format a translation warning as a single tab-separated line."""
    fields = [warning.type.value, warning.entity_id or "", warning.message]
    return "\t".join(" ".join(f.split()) for f in fields)


@app.command()
def main(
    input_file: Annotated[
        Path,
        typer.Argument(
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            help="Minerva model JSON file.",
        ),
    ],
    warnings_only: Annotated[
        bool,
        typer.Option(
            help="If set, only write translation warnings; do not write the GO-CAM model.",
        ),
    ] = False,
):
    """
    Translate a Minerva model JSON file to a GO-CAM model and report translation warnings.

    Development aid only; this is not a pipeline step.
    """
    with open(input_file, "r") as f:
        minerva_model = json.load(f)

    translation_result = MinervaWrapper.translate(minerva_model)

    if not warnings_only:
        print(translation_result.result.model_dump_json(indent=2, exclude_none=True))

    for warning in sorted(
        translation_result.warnings,
        key=lambda w: (w.type.value, w.entity_id or "", w.message),
    ):
        print(format_warning(warning), file=sys.stderr)


if __name__ == "__main__":
    app()
