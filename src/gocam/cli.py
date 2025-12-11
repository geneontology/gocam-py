import csv
import datetime
import io
import json
import logging
import os
import sys
import tarfile
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from enum import Enum
from pathlib import Path
from typing import List, Optional
from urllib.request import urlretrieve

import typer
import yaml
from typing_extensions import Annotated, Literal

from gocam import __version__
from gocam.datamodel import Model
from gocam.indexing.Flattener import Flattener
from gocam.indexing.Indexer import Indexer
from gocam.translation import MinervaWrapper
from gocam.translation.cx2 import model_to_cx2
from gocam.translation.networkx.model_network_translator import ModelNetworkTranslator


logger = logging.getLogger(__name__)


class TranslationFormat(str, Enum):
    """Translation format options."""

    NETWORKX = "networkx"
    CX2 = "cx2"


app = typer.Typer(
    help="A CLI for interacting with GO-CAMs.", rich_markup_mode="markdown"
)


def version_callback(value: bool):
    if value:
        typer.echo(__version__)
        raise typer.Exit()


@app.callback()
def main(
    verbose: Annotated[
        int, typer.Option("--verbose", "-v", count=True, help="Increase verbosity")
    ] = 0,
    quiet: Annotated[
        bool, typer.Option("--quiet", "-q", help="Suppress output")
    ] = False,
    stacktrace: Annotated[
        bool,
        typer.Option(help="If set then show full stacktrace on error"),
    ] = False,
    version: Annotated[
        Optional[bool],
        typer.Option(
            callback=version_callback,
            is_eager=True,
            help="Show version and exit",
        ),
    ] = None,
):
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


@app.command()
def fetch(
    model_ids: Annotated[
        Optional[List[str]], typer.Argument(help="Model IDs to fetch")
    ] = None,
    format: Annotated[
        Literal["json", "yaml"], typer.Option("--format", "-f", help="Input format")
    ] = "yaml",
    add_indexes: Annotated[
        bool,
        typer.Option(
            help="Add indexes (closures, counts) to the model",
        ),
    ] = False,
    as_minerva: Annotated[
        bool,
        typer.Option(help="Export as minerva json/yaml"),
    ] = False,
):
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
            typer.echo(json.dumps(model_dict, indent=2))
        elif format == "yaml":
            typer.echo("---")
            typer.echo(yaml.dump(model_dict, sort_keys=False))
        else:
            typer.echo(model_dict)


@app.command()
def convert(
    output_format: Annotated[
        Literal["cx2", "owl"],
        typer.Option("--output-format", "-O", help="Output format"),
    ],
    model: Annotated[typer.FileText, typer.Argument(help="Input model file")] = "-",
    input_format: Annotated[
        Optional[Literal["json", "yaml"]],
        typer.Option(
            "--input-format",
            "-I",
            help="Input format. Not required unless reading from stdin.",
        ),
    ] = None,
    output: Annotated[
        typer.FileTextWrite, typer.Option("--output", "-o", help="Output file")
    ] = "-",
    dot_layout: Annotated[
        bool, typer.Option(help="Apply dot layout (requires Graphviz)")
    ] = False,
    ndex_upload: Annotated[
        bool, typer.Option(help="Upload to NDEx (only for CX2)")
    ] = False,
):
    """Convert GO-CAM models.

    Currently supports converting to CX2 format and uploading to NDEx.
    """
    if ndex_upload and output_format != "cx2":
        raise typer.BadParameter("NDEx upload requires output format to be CX2")

    if input_format is None:
        if model.name.endswith(".json"):
            input_format = "json"
        elif model.name.endswith(".yaml"):
            input_format = "yaml"
        else:
            raise typer.BadParameter("Could not infer input format")

    if input_format == "json":
        # Use Pydantic's built-in JSON parser for better performance
        json_content = model.read()
        try:
            # Try to parse as a single model first
            models = [Model.model_validate_json(json_content)]
        except Exception:
            # If that fails, try parsing as a list of models
            deserialized = json.loads(json_content)
            if isinstance(deserialized, list):
                models = [Model.model_validate(m) for m in deserialized]
            else:
                raise
    elif input_format == "yaml":
        deserialized = list(yaml.safe_load_all(model))
        models = [Model.model_validate(m) for m in deserialized]
    else:
        raise typer.BadParameter("Invalid input format")

    try:
        logger.debug(f"Parsed {len(models)} models from input")
    except Exception as e:
        raise typer.BadParameter(f"Could not load model: {e}")

    if output_format == "cx2":
        if len(models) != 1:
            raise typer.BadParameter("CX2 format only supports a single model")
        model_obj = models[0]
        cx2 = model_to_cx2(model_obj, apply_dot_layout=dot_layout)

        if ndex_upload:
            try:
                import ndex2
            except ImportError:
                raise typer.BadParameter(
                    "ndex2 package is required for NDEx upload. Install with: pip install ndex2"
                )

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

            typer.echo(
                f"View network at: 'https://www.ndexbio.org/viewer/networks/{network_id}"
            )
        else:
            typer.echo(json.dumps(cx2), file=output)
    elif output_format == "owl":
        from gocam.translation.tbox_translator import TBoxTranslator

        tbox_translator = TBoxTranslator()
        tbox_translator.load_models(models)
        tbox_translator.save_ontology(output.name, serialization="ofn")


@app.command()
def index_models(
    input_file: Annotated[
        Path, typer.Argument(help="Input file containing models", exists=True)
    ],
    input_format: Annotated[
        Optional[Literal["json", "yaml"]],
        typer.Option(
            "--input-format",
            "-I",
            help="Input format. Not required unless reading from stdin or file has no extension.",
        ),
    ] = None,
    output_format: Annotated[
        Literal["json", "yaml"],
        typer.Option(
            "--output-format", "-O", help="Output format for the indexed models."
        ),
    ] = "yaml",
    output_file: Annotated[
        Optional[Path],
        typer.Option(
            "--output-file",
            "-o",
            help="Output file. If not specified, write to stdout.",
        ),
    ] = None,
    reindex: Annotated[
        bool,
        typer.Option(help="Reindex models that already have indexes"),
    ] = False,
):
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
            raise typer.BadParameter(
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
                raise typer.BadParameter("JSON input must be a list of models")
            for model_dict in data:
                try:
                    # Use Pydantic's built-in JSON parser for better performance
                    model_json = json.dumps(model_dict)
                    model = Model.model_validate_json(model_json)
                    models.append(model)
                except Exception as e:
                    typer.echo(f"Warning: Could not load model: {e}", err=True)
    else:  # yaml
        # For YAML, support multiple documents
        with open(input_path, "r") as f:
            yaml_content = f.read()
        for doc in yaml.safe_load_all(yaml_content):
            try:
                model = Model.model_validate(doc)
                models.append(model)
            except Exception as e:
                typer.echo(f"Warning: Could not load model: {e}", err=True)

    typer.echo(f"Loaded {len(models)} models from {input_file}", err=True)

    # Index models
    indexer = Indexer()
    for model in models:
        try:
            indexer.index_model(model, reindex=reindex)
        except Exception as e:
            typer.echo(f"Warning: Could not index model {model.id}: {e}", err=True)

    typer.echo(f"Indexed {len(models)} models", err=True)

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
        typer.echo(f"Wrote indexed models to {output_file}", err=True)
    else:
        typer.echo(output_content)


@app.command()
def flatten_models(
    input_file: Annotated[
        Path, typer.Argument(help="Input file containing indexed models", exists=True)
    ],
    input_format: Annotated[
        Optional[Literal["json", "yaml"]],
        typer.Option(
            "--input-format",
            "-I",
            help="Input format. Not required unless reading from stdin or file has no extension.",
        ),
    ] = None,
    output_format: Annotated[
        Literal["json", "jsonl", "tsv"],
        typer.Option(
            "--output-format", "-O", help="Output format for the flattened data."
        ),
    ] = "jsonl",
    output_file: Annotated[
        Optional[Path],
        typer.Option(
            "--output-file",
            "-o",
            help="Output file. If not specified, write to stdout.",
        ),
    ] = None,
    fields: Annotated[
        Optional[str],
        typer.Option(
            "--fields",
            "-f",
            help="Comma-separated list of fields to include in the output. If not specified, all fields are included.",
        ),
    ] = None,
):
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
            raise typer.BadParameter(
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
                raise typer.BadParameter("JSON input must be a list of models")
            for model_dict in data:
                try:
                    # Use Pydantic's built-in JSON parser for better performance
                    model_json = json.dumps(model_dict)
                    model = Model.model_validate_json(model_json)
                    models.append(model)
                except Exception as e:
                    typer.echo(f"Warning: Could not load model: {e}", err=True)
    else:  # yaml
        # For YAML, support multiple documents
        with open(input_path, "r") as f:
            yaml_content = f.read()
        for doc in yaml.safe_load_all(yaml_content):
            try:
                model = Model.model_validate(doc)
                models.append(model)
            except Exception as e:
                typer.echo(f"Warning: Could not load model: {e}", err=True)

    typer.echo(f"Loaded {len(models)} models from {input_file}", err=True)

    # Flatten models
    flattener = Flattener(fields=field_list)
    rows = []
    for model in models:
        try:
            row = flattener.flatten(model)
            rows.append(row)
        except Exception as e:
            typer.echo(f"Warning: Could not flatten model {model.id}: {e}", err=True)

    typer.echo(f"Flattened {len(rows)} models", err=True)

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
            output_buffer = io.StringIO()
            writer = csv.DictWriter(
                output_buffer,
                fieldnames=fieldnames,
                delimiter="\t",
                extrasaction="ignore",
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
        typer.echo(f"Wrote flattened models to {output_file}", err=True)
    else:
        typer.echo(output_content)


@app.command("translate-collection")
def translate_collection(
    url: Annotated[
        str, typer.Option(help="URL to download the GO-CAM models tarball")
    ] = "http://current.geneontology.org/products/json/noctua-models-json.tgz",
    format: Annotated[
        Optional[List[TranslationFormat]],
        typer.Option(
            "--format", help="Formats to translate models to (can specify multiple)"
        ),
    ] = None,
    output: Annotated[
        str,
        typer.Option(
            help="Base output directory (format subdirectories will be created automatically)",
        ),
    ] = "/tmp/gocam-output",
    limit: Annotated[
        Optional[int],
        typer.Option(help="Limit the number of models to process (for testing)"),
    ] = None,
    archive: Annotated[
        bool,
        typer.Option(
            help="Create gzipped tar archives of output directories",
        ),
    ] = True,
    max_workers: Annotated[
        int,
        typer.Option(
            help="Maximum number of concurrent workers for parallel processing",
        ),
    ] = 10,
    batch_size: Annotated[
        int,
        typer.Option(
            help="Number of models to process in each batch (0 for no batching)",
        ),
    ] = 100,
):
    """
    Downloads a tarball of GO-CAM models in minerva JSON format, then translates
    each model through networkx and/or cx2 formats.

    Examples:

        # Translate all models to both formats
        gocam translate-collection

        # Only NetworkX format with custom output directory
        gocam translate-collection --format networkx --output ./my_output

        # Test with limited models
        gocam translate-collection --limit 5
    """

    def download_and_extract_tarball(url: str, extract_dir: str) -> str:
        """Download and extract a gzipped tarball."""
        logger.info(f"Downloading tarball from {url}")

        with tempfile.NamedTemporaryFile(suffix=".tgz") as tmp_file:
            urlretrieve(url, tmp_file.name)

            logger.info(f"Extracting tarball to {extract_dir}")
            with tarfile.open(tmp_file.name, "r:gz") as tar:
                tar.extractall(path=extract_dir)

                # Find the extracted directory
                extracted_dirs = [
                    d
                    for d in os.listdir(extract_dir)
                    if os.path.isdir(os.path.join(extract_dir, d))
                ]
                if extracted_dirs:
                    return os.path.join(extract_dir, extracted_dirs[0])
                else:
                    return extract_dir

    def find_json_models(directory: str, limit: Optional[int] = None) -> List[str]:
        """Find JSON model files in a directory, optionally limiting the number."""
        json_files = []
        for root, dirs, files in os.walk(directory):
            for file in files:
                if file.endswith(".json"):
                    json_files.append(os.path.join(root, file))
                    if limit and len(json_files) >= limit:
                        logger.info(
                            f"Found {len(json_files)} JSON model files (limited)"
                        )
                        return json_files

        logger.info(f"Found {len(json_files)} JSON model files")
        return json_files

    def translate_to_networkx(model: Model, output_path: str, model_filename: str):
        """Translate a model to NetworkX format and save as JSON."""
        try:
            # No indexer necessary - causal associations are already populated by minerva_wrapper
            translator = ModelNetworkTranslator()
            g2g_graph = translator.translate_models([model])

            # Check if the graph has any edges before writing
            if g2g_graph.number_of_edges() > 0:
                # Convert graph to dict directly to avoid redundant translation
                g2g_dict = translator._graph_to_dict(
                    g2g_graph, [model], include_model_info=True
                )

                # Save to file
                output_file = os.path.join(
                    output_path, f"{model_filename}_networkx.json"
                )
                with open(output_file, "w") as f:
                    json.dump(g2g_dict, f)

        except Exception as e:
            logger.error(f"Failed to translate {model_filename} to NetworkX: {e}")

    def translate_to_cx2(model: Model, output_path: str, model_filename: str):
        """Translate a model to CX2 format and save as JSON."""
        try:
            # Check if model has causal associations before translating
            has_causal_edges = False
            if model.activities:
                for activity in model.activities:
                    if activity.causal_associations:
                        has_causal_edges = True
                        break

            # Only create file if there are causal edges (matching NetworkX behavior)
            if has_causal_edges:
                # Translate the model
                cx2_data = model_to_cx2(
                    model, validate_iquery_gene_symbol_pattern=False
                )
                # Save to file
                output_file = os.path.join(output_path, f"{model_filename}_cx2.json")
                with open(output_file, "w") as f:
                    json.dump(cx2_data, f)
        except Exception as e:
            logger.error(f"Failed to translate {model_filename} to CX2: {e}")

    # Set default formats if not specified
    if format is None:
        format = [TranslationFormat.NETWORKX, TranslationFormat.CX2]

    # Ensure format is a list (Click's multiple option sometimes returns a tuple)
    if isinstance(format, tuple):
        format = list(format)

    # Create output directories based on selected formats
    output_paths = {}
    for fmt in format:
        # Always create subdirectories for each format
        output_paths[fmt.value] = os.path.join(output, fmt.value)
        os.makedirs(output_paths[fmt.value], exist_ok=True)
        typer.echo(
            f"Created output directory for {fmt.value}: {output_paths[fmt.value]}",
            err=True,
        )

    # Create temporary directory for extraction
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Download and extract the tarball
            extracted_dir = download_and_extract_tarball(url, temp_dir)

            # Find all JSON model files
            json_files = find_json_models(extracted_dir, limit)

            # Progress tracking
            progress_lock = threading.Lock()
            processed_count = 0
            failed_count = 0
            processed_model_ids = []

            def report_progress(phase, count, total):
                percentage = (count / total) * 100
                typer.echo(f"{phase}: {count}/{total} ({percentage:.1f}%)", err=True)

            typer.echo(
                f"Processing {len(json_files)} models (load+translate)...", err=True
            )

            def load_and_translate_single(file_path: str) -> bool:
                """Load a model from file and translate it to all requested formats."""
                nonlocal processed_count, failed_count, processed_model_ids

                model_filename = os.path.splitext(os.path.basename(file_path))[0]
                try:
                    # Load model
                    with open(file_path, "r") as f:
                        minerva_dict = json.load(f)
                        model = MinervaWrapper.minerva_object_to_model(minerva_dict)

                    # Translate to all requested formats
                    for format_name in format:
                        if format_name == TranslationFormat.NETWORKX:
                            translate_to_networkx(
                                model, output_paths["networkx"], model_filename
                            )
                        elif format_name == TranslationFormat.CX2:
                            translate_to_cx2(model, output_paths["cx2"], model_filename)

                    # Track success
                    with progress_lock:
                        processed_count += 1
                        processed_model_ids.append(model.id)
                        if (
                            processed_count % 100 == 0
                            or processed_count % max(1, len(json_files) // 20) == 0
                        ):
                            report_progress(
                                "Processed", processed_count, len(json_files)
                            )

                    return True

                except Exception as e:
                    logger.error(f"Failed to process model from {file_path}: {e}")
                    with progress_lock:
                        failed_count += 1
                    return False

            # Process all files concurrently
            stop_event = threading.Event()

            def interruptible_translate(json_file):
                if stop_event.is_set():
                    return False
                return load_and_translate_single(json_file)

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [
                    executor.submit(interruptible_translate, json_file)
                    for json_file in json_files
                ]

                try:
                    for future in as_completed(futures):
                        if stop_event.is_set():
                            break
                        future.result()  # This will raise any exceptions
                except KeyboardInterrupt:
                    logger.info("Interrupted by user, cancelling remaining tasks...")
                    stop_event.set()
                    executor.shutdown(wait=False, cancel_futures=True)
                    raise

            typer.echo(
                f"Processing complete. Successfully processed: {processed_count}, Failed: {failed_count}",
                err=True,
            )

            # Create tar.gz archives if requested
            if archive:

                def create_archive(directory_path: str, format_name: str):
                    """Create a gzipped tar archive of a directory."""
                    if not os.path.exists(directory_path) or not os.listdir(
                        directory_path
                    ):
                        typer.echo(
                            f"Skipping archive for {format_name}: directory is empty or doesn't exist",
                            err=True,
                        )
                        return

                    # Create archive filename
                    schema_version = __version__
                    archive_name = f"go-cam-{format_name}.tar.gz"
                    # Place archive in the base output directory
                    archive_path = os.path.join(output, archive_name)

                    typer.echo(
                        f"Creating {format_name} archive: {archive_path}", err=True
                    )

                    try:
                        # Create metadata file
                        metadata = {
                            "schema_version": schema_version,
                            "format": format_name,
                            "created_at": datetime.datetime.now().isoformat(),
                            "source_url": url,
                            "models_processed": processed_count,
                            "models_failed": failed_count,
                        }

                        metadata_path = os.path.join(
                            directory_path, "gocam_metadata.json"
                        )
                        with open(metadata_path, "w") as f:
                            json.dump(metadata, f, indent=2)

                        # Create model index file
                        model_index = {
                            "model_ids": processed_model_ids,
                            "total_models": len(processed_model_ids),
                        }

                        index_path = os.path.join(
                            directory_path, "gocam_model_index.json"
                        )
                        with open(index_path, "w") as f:
                            json.dump(model_index, f, indent=2)

                        with tarfile.open(archive_path, "w:gz") as tar:
                            # Add all files in the directory, including metadata
                            for file_name in os.listdir(directory_path):
                                file_path = os.path.join(directory_path, file_name)
                                if os.path.isfile(file_path):
                                    tar.add(file_path, arcname=file_name)

                        # Clean up temporary files
                        os.unlink(metadata_path)
                        os.unlink(index_path)

                        typer.echo(
                            f"Successfully created archive: {archive_path}", err=True
                        )

                    except Exception as e:
                        logger.error(f"Failed to create archive for {format_name}: {e}")

                # Create archives for each format that was processed
                if "networkx" in format and processed_count > 0:
                    create_archive(output_paths["networkx"], "networkx")

                if "cx2" in format and processed_count > 0:
                    create_archive(output_paths["cx2"], "cx2")

        except Exception as e:
            logger.error(f"Translation failed with error: {e}")
            raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
