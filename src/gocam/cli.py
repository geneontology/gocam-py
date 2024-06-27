import json
import logging
import sys
import warnings
from typing import Optional

import click
import yaml
from linkml_runtime.dumpers import json_dumper
from pydantic import BaseModel

from gocam.translation.minerva_wrapper import MinervaWrapper

index_type_option = click.option(
    "--index-type",
    "-t",
    default="simple",
    show_default=True,
    help="Type of index to create. Values: simple, llm",
)

logger = logging.getLogger(__name__)

warnings.filterwarnings("ignore", module="duckdb_engine")

def get_wrapper():
    return MinervaWrapper()


format_choice = click.Choice(["json", "yaml", "tsv"])


include_internal_option = click.option("--include-internal/--no-include-internal", default=False, show_default=True)


@click.group()
@click.option("-v", "--verbose", count=True)
@click.option("-q", "--quiet/--no-quiet")
@click.option(
    "--stacktrace/--no-stacktrace",
    default=False,
    show_default=True,
    help="If set then show full stacktrace on error",
)
@click.pass_context
def cli(ctx, verbose: int, quiet: bool, stacktrace: bool):
    """A CLI for interacting with the linkml-store."""
    if not stacktrace:
        sys.tracebacklimit = 0
    logger = logging.getLogger()
    # Set handler for the root logger to output to the console
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))

    # Clear existing handlers to avoid duplicate messages if function runs multiple times
    logger.handlers = []

    # Add the newly created console handler to the logger
    logger.addHandler(console_handler)
    if verbose >= 2:
        logger.setLevel(logging.DEBUG)
    elif verbose == 1:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
    if quiet:
        logger.setLevel(logging.ERROR)


@cli.command()
@click.option("--format", "-f", type=format_choice, help="Input format")
@click.argument("model_ids", nargs=-1)
def fetch(model_ids, format):
    wrapper = get_wrapper()
    if not model_ids:
        model_ids = wrapper.models_ids()
    for model_id in model_ids:
        model = wrapper.object_by_id(model_id)
        model_dict = model.model_dump(exclude_none=True, exclude_unset=True)
        if format is None:
            format = "yaml"
        if format == "json":
            print(json.dumps(model_dict, indent=2))
        elif format == "yaml":
            print("---")
            print(yaml.dump(model_dict, sort_keys=False))
        else:
            click.echo(model.model_dump())


if __name__ == "__main__":
    cli()
