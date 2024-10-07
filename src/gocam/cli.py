import json
import logging
import os
import sys

import click
import ndex2
import yaml

from gocam import __version__
from gocam.datamodel import Model
from gocam.translation import MinervaWrapper
from gocam.translation.cx2 import model_to_cx2


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
        model_dict = model.model_dump(exclude_none=True)

        if format == "json":
            click.echo(json.dumps(model_dict, indent=2))
        elif format == "yaml":
            click.echo("---")
            click.echo(yaml.dump(model_dict, sort_keys=False))
        else:
            click.echo(model.model_dump())


@cli.command()
@click.option(
    "--input-format",
    "-I",
    type=click.Choice(["json", "yaml"]),
    help="Input format. Not required unless reading from stdin.",
)
@click.option("--output-format", "-O", type=click.Choice(["cx2"]), required=True)
@click.option("--output", "-o", type=click.File("w"), default="-")
@click.option("--ndex-upload", is_flag=True, help="Upload to NDEx")
@click.argument("model", type=click.File("r"), default="-")
def convert(model, input_format, output_format, output, ndex_upload):
    """Convert GO-CAM models."""
    if ndex_upload and output_format != "cx2":
        raise click.UsageError("NDEx upload requires output format to be CX2")

    if input_format is None:
        if model.name.endswith(".json"):
            input_format = "json"
        elif model.name.endswith(".yaml"):
            input_format = "yaml"
        else:
            raise click.BadParameter("Could not infer input format")

    if input_format == "json":
        deserialized = json.load(model)
    elif input_format == "yaml":
        deserialized = yaml.safe_load(model)
    else:
        raise click.UsageError("Invalid input format")

    try:
        model = Model.model_validate(deserialized)
    except Exception as e:
        raise click.UsageError(f"Could not load model: {e}")

    if output_format == "cx2":
        cx2 = model_to_cx2(model)

        if ndex_upload:
            # This is very basic proof-of-concept usage of the NDEx client. Once we have a better
            # idea of how we want to use it, we can refactor this to allow more CLI options for
            # connection details, visibility, adding the new network to a group, etc. At that point
            # we can also consider moving upload functionality to a separate command.
            client = ndex2.client.Ndex2(
                host=os.getenv("NDEX_HOST"),
                username=os.getenv("NDEX_USERNAME"),
                password=os.getenv("NDEX_PASSWORD"),
            )
            url = client.save_new_cx2_network(cx2, visibility="PRIVATE")
            network_id = url.rsplit("/", 1)[-1]

            # Make the network searchable
            client.set_network_system_properties(network_id, {"index_level": "META"})

            click.echo(
                f"View network at: 'https://www.ndexbio.org/viewer/networks/{network_id}"
            )
        else:
            click.echo(json.dumps(cx2), file=output)


if __name__ == "__main__":
    cli()
