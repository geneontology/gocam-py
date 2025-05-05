import json
import logging
import os
import sys
from pathlib import Path
from typing import List, Optional

import click
import ndex2
import yaml
from linkml_runtime.loaders import json_loader, yaml_loader

from gocam import __version__
from gocam.datamodel import Model
from gocam.translation import MinervaWrapper
from gocam.translation.cx2 import model_to_cx2
from gocam.indexing.indexer import Indexer


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
@click.option(
    "--add-indexes/--no-add-indexes",
    default=False,
    show_default=True,
    help="Add indexes (closures, counts) to the model",
)
@click.argument("model_ids", nargs=-1)
def fetch(model_ids, format, add_indexes):
    """Fetch GO-CAM models."""
    wrapper = MinervaWrapper()
    indexer = None
    if add_indexes:
        indexer = Indexer()

    if not model_ids:
        model_ids = wrapper.models_ids()

    for model_id in model_ids:
        model = wrapper.fetch_model(model_id)
        if indexer:
            indexer.index_model(model)
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
@click.option("--dot-layout", is_flag=True, help="Apply dot layout (requires Graphviz)")
@click.option("--ndex-upload", is_flag=True, help="Upload to NDEx (only for CX2)")
@click.argument("model", type=click.File("r"), default="-")
def convert(model, input_format, output_format, output, dot_layout, ndex_upload):
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
        cx2 = model_to_cx2(model, apply_dot_layout=dot_layout)

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


@cli.command()
@click.option(
    "--input-format",
    "-I",
    type=click.Choice(["json", "yaml"]),
    help="Input format. Not required unless reading from stdin or file has no extension.",
)
@click.option(
    "--output-format",
    "-O",
    type=click.Choice(["json", "yaml"]),
    default="yaml",
    help="Output format for the indexed models.",
)
@click.option(
    "--output-file",
    "-o",
    type=click.Path(writable=True),
    help="Output file. If not specified, write to stdout.",
)
@click.option(
    "--reindex/--no-reindex",
    default=False,
    show_default=True,
    help="Reindex models that already have indexes",
)
@click.argument("input_file", type=click.Path(exists=True))
def index_models(input_file, input_format, output_format, output_file, reindex):
    """
    Index a collection of GO-CAM models.

    This command takes a file containing a list of GO-CAM models (in JSON or YAML format),
    adds indexes to each model, and outputs the indexed models.

    For YAML input, the file can contain multiple documents separated by '---'.
    """
    input_path = Path(input_file)

    # Determine input format if not specified
    if input_format is None:
        if input_path.suffix.lower() == ".json":
            input_format = "json"
        elif input_path.suffix.lower() in [".yaml", ".yml"]:
            input_format = "yaml"
        else:
            raise click.BadParameter(
                "Could not infer input format from file extension. Please specify --input-format."
            )

    # Load models
    models: List[Model] = []
    if input_format == "json":
        # For JSON, expect a list of model objects
        with open(input_path, "r") as f:
            data = json.load(f)
            if not isinstance(data, list):
                raise click.BadParameter("JSON input must be a list of models")
            for model_dict in data:
                try:
                    model = Model.model_validate(model_dict)
                    models.append(model)
                except Exception as e:
                    click.echo(f"Warning: Could not load model: {e}", err=True)
    else:  # yaml
        # For YAML, support multiple documents
        with open(input_path, "r") as f:
            yaml_content = f.read()
        for doc in yaml.safe_load_all(yaml_content):
            try:
                model = Model.model_validate(doc)
                models.append(model)
            except Exception as e:
                click.echo(f"Warning: Could not load model: {e}", err=True)

    click.echo(f"Loaded {len(models)} models from {input_file}", err=True)

    # Index models
    indexer = Indexer()
    for model in models:
        try:
            indexer.index_model(model, reindex=reindex)
        except Exception as e:
            click.echo(f"Warning: Could not index model {model.id}: {e}", err=True)

    click.echo(f"Indexed {len(models)} models", err=True)

    # Output indexed models
    if output_format == "json":
        output_data = [model.model_dump(exclude_none=True) for model in models]
        output_content = json.dumps(output_data, indent=2)
    else:  # yaml
        # For YAML, output multiple documents
        output_content = ""
        for model in models:
            output_content += "---\n"
            output_content += yaml.dump(
                model.model_dump(exclude_none=True), sort_keys=False
            )

    if output_file:
        with open(output_file, "w") as f:
            f.write(output_content)
        click.echo(f"Wrote indexed models to {output_file}", err=True)
    else:
        click.echo(output_content)


if __name__ == "__main__":
    cli()
