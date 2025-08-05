import csv
import json
import logging
import os
import sys
from pathlib import Path
from typing import List, Optional

import click
import yaml
from linkml_runtime.loaders import json_loader, yaml_loader
from mkdocs.commands.serve import serve

from gocam import __version__
from gocam.datamodel import Model
from gocam.translation import MinervaWrapper
from gocam.translation.cx2 import model_to_cx2
from gocam.indexing.Indexer import Indexer
from gocam.indexing.Flattener import Flattener

import logging

logger = logging.getLogger(__name__)

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
@click.option(
    "--as-minerva/--no-as-minerva",
    default=False,
    show_default=True,
    help="Export as minerva json/yaml",
)
@click.argument("model_ids", nargs=-1)

def fetch(model_ids, as_minerva, format, add_indexes):
    """Fetch GO-CAM models.

    TODO: this currently fetches from a pre-filtered set of models.

    Fetch and convert to GO-CAM yaml:

        gocam fetch 61e0e55600000624

    Add indexes (closures, counts) to the model:

        gocam fetch --add-indexes 61e0e55600000624

    Fetch, preserving minerva low-level format:

        gocam fetch --as-minerva gomodel:YeastPathways_LYSDEGII-PWY

    Note: this should be used mostly for debugging purposes.

    """
    wrapper = MinervaWrapper()
    indexer = None
    if add_indexes:
        indexer = Indexer()

    if not model_ids:
        model_ids = wrapper.models_ids()

    for model_id in model_ids:
        if as_minerva:
            model_dict = wrapper.fetch_minerva_object(model_id)
        else:
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
            click.echo(model_dict)


@cli.command()
@click.option(
    "--input-format",
    "-I",
    type=click.Choice(["json", "yaml"]),
    help="Input format. Not required unless reading from stdin.",
)
@click.option("--output-format", "-O", type=click.Choice(["cx2", "owl"]), required=True)
@click.option("--output", "-o", type=click.File("w"), default="-")
@click.option("--dot-layout", is_flag=True, help="Apply dot layout (requires Graphviz)")
@click.option("--ndex-upload", is_flag=True, help="Upload to NDEx (only for CX2)")
@click.argument("model", type=click.File("r"), default="-")
def convert(model, input_format, output_format, output, dot_layout, ndex_upload):
    """Convert GO-CAM models.

    Currently supports converting to CX2 format and uploading to NDEx.
    """
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
        # Use Pydantic's built-in JSON parser for better performance
        json_content = model.read()
        try:
            # Try to parse as a single model first
            models = [Model.model_validate_json(json_content)]
            logger.info("Parsing a single model from input")
        except Exception:
            # If that fails, try parsing as a list of models
            deserialized = json.loads(json_content)
            if isinstance(deserialized, list):
                logger.info(f"Parsing {len(deserialized)} models from input")
                models = [Model.model_validate(m) for m in deserialized]
            else:
                raise
    elif input_format == "yaml":
        deserialized = list(yaml.safe_load_all(model))
        logger.info(f"Parsing {len(deserialized)} models from input")
        models = [Model.model_validate(m) for m in deserialized]
    else:
        raise click.UsageError("Invalid input format")

    try:
        logger.info(f"Parsed {len(models)} models from input")
    except Exception as e:
        raise click.UsageError(f"Could not load model: {e}")

    if output_format == "cx2":
        if len(models) != 1:
            raise click.UsageError("CX2 format only supports a single model")
        model = models[0]
        cx2 = model_to_cx2(model, apply_dot_layout=dot_layout)

        if ndex_upload:
            import ndex2

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
    elif output_format == "owl":
        from gocam.translation.tbox_translator import TBoxTranslator
        tbox_translator = TBoxTranslator()
        tbox_translator.load_models(models)
        tbox_translator.save_ontology(output.name, serialization="ofn")


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
        with open(input_path, "rb") as f:
            json_content = f.read()
            data = json.loads(json_content)
            if not isinstance(data, list):
                raise click.BadParameter("JSON input must be a list of models")
            for model_dict in data:
                try:
                    # Use Pydantic's built-in JSON parser for better performance
                    model_json = json.dumps(model_dict)
                    model = Model.model_validate_json(model_json)
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
    type=click.Choice(["json", "jsonl", "tsv"]),
    default="jsonl",
    help="Output format for the flattened data.",
)
@click.option(
    "--output-file",
    "-o",
    type=click.Path(writable=True),
    help="Output file. If not specified, write to stdout.",
)
@click.option(
    "--fields",
    "-f",
    help="Comma-separated list of fields to include in the output. If not specified, all fields are included.",
)
@click.argument("input_file", type=click.Path(exists=True))
def flatten_models(input_file, input_format, output_format, output_file, fields):
    """
    Flatten indexed GO-CAM models into tabular format.

    This command takes a file containing indexed GO-CAM models (in JSON or YAML format),
    flattens them into rows, and outputs them in JSON, JSONL, or TSV format.

    The models should already be indexed (have query_index populated).
    Use the index-models command first if needed.

    Example:
        gocam flatten-models -O tsv indexed-models.yaml -o models.tsv
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

    # Parse fields if provided
    field_list = None
    if fields:
        field_list = [f.strip() for f in fields.split(",")]

    # Load models
    models: List[Model] = []
    if input_format == "json":
        # For JSON, expect a list of model objects
        with open(input_path, "rb") as f:
            json_content = f.read()
            data = json.loads(json_content)
            if not isinstance(data, list):
                raise click.BadParameter("JSON input must be a list of models")
            for model_dict in data:
                try:
                    # Use Pydantic's built-in JSON parser for better performance
                    model_json = json.dumps(model_dict)
                    model = Model.model_validate_json(model_json)
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

    # Flatten models
    flattener = Flattener(fields=field_list)
    rows = []
    for model in models:
        try:
            row = flattener.flatten(model)
            rows.append(row)
        except Exception as e:
            click.echo(f"Warning: Could not flatten model {model.id}: {e}", err=True)

    click.echo(f"Flattened {len(rows)} models", err=True)

    # Prepare output
    output_content = None
    
    if output_format == "json":
        output_content = json.dumps(rows, indent=2)
    elif output_format == "jsonl":
        output_lines = [json.dumps(row) for row in rows]
        output_content = "\n".join(output_lines)
    elif output_format == "tsv":
        if not rows:
            output_content = ""
        else:
            # Get all unique field names
            all_fields = set()
            for row in rows:
                all_fields.update(row.keys())
            fieldnames = sorted(all_fields)
            
            # Convert to TSV
            import io
            output_buffer = io.StringIO()
            writer = csv.DictWriter(
                output_buffer, 
                fieldnames=fieldnames, 
                delimiter="\t",
                extrasaction='ignore'
            )
            writer.writeheader()
            for row in rows:
                # Convert list values to comma-separated strings
                tsv_row = {}
                for k, v in row.items():
                    if isinstance(v, list):
                        tsv_row[k] = ",".join(str(item) for item in v)
                    else:
                        tsv_row[k] = v
                writer.writerow(tsv_row)
            output_content = output_buffer.getvalue()

    # Write output
    if output_file:
        with open(output_file, "w") as f:
            f.write(output_content)
        click.echo(f"Wrote flattened models to {output_file}", err=True)
    else:
        click.echo(output_content)


if __name__ == "__main__":
    cli()
