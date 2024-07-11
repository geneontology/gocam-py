import json
import logging
import sys

import click
import yaml

from gocam import __version__
from gocam.translation import MinervaWrapper


@click.group()
@click.option("-v", "--verbose", count=True)
@click.option("-q", "--quiet/--no-quiet")
@click.option(
    "--stacktrace/--no-stacktrace",
    default=False,
    show_default=True,
    help="If set then show full stacktrace on error",
)
@click.version_option(__version__)
def cli(verbose: int, quiet: bool, stacktrace: bool):
    """A CLI for interacting with GO-CAMs."""
    if not stacktrace:
        sys.tracebacklimit = 0

    logger = logging.getLogger()
    # Set handler for the root logger to output to the console
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(
        logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    )

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
@click.option(
    "--format",
    "-f",
    type=click.Choice(["json", "yaml"]),
    default="yaml",
    show_default=True,
    help="Input format",
)
@click.argument("model_ids", nargs=-1)
def fetch(model_ids, format):
    """Fetch GO-CAM models."""
    wrapper = MinervaWrapper()

    if not model_ids:
        model_ids = wrapper.models_ids()

    for model_id in model_ids:
        model = wrapper.fetch_model(model_id)
        model_dict = model.model_dump(exclude_none=True, exclude_unset=True)

        if format == "json":
            click.echo(json.dumps(model_dict, indent=2))
        elif format == "yaml":
            click.echo("---")
            click.echo(yaml.dump(model_dict, sort_keys=False))
        else:
            click.echo(model.model_dump())


if __name__ == "__main__":
    cli()
