# Pipeline Scripts

This directory contains scripts for processing steps in the GOC pipeline.

## Setup

These scripts are not included in the published `gocam-py` package, so you will need to set up a
local development environment to run them. Use `uv` to install the necessary dependencies:

```bash
uv sync --all-extras
```

## convert_minerva_models_to_gocam_models.py

Converts Minerva model JSON files to GO-CAM model format, doing basic filtering to eliminate obvious
non-GO-CAM models. The output files are intended for internal QC analysis and further filtering.

**Inputs:**

- `--input-dir`: Directory containing Minerva model JSON files
- `--output-dir`: Directory to save converted GO-CAM models (required unless using `--dry-run`)

**Filtering:**

Models are filtered out and not written if they:

- Have no activity edges (no activities linked by causal associations or shared chemical
  inputs/outputs)
- Use complement types

**Outputs:**

- GO-CAM model JSON files (one per input file) in the specified output directory

**Options:**

- `--dry-run`: Perform conversion without writing output files
- `--verbose` / `-v`: Increase verbosity (use multiple times for more detail: `-v` for INFO, `-vv`
  for DEBUG)
- `--limit N`: Limit processing to first N models (0 = no limit)

**Usage:**

```bash
python pipeline/convert_minerva_models_to_gocam_models.py --input-dir /path/to/minerva/models --output-dir /path/to/gocam/models
```

## filter_true_gocam_models.py

Filters a collection of models based on whether they meet the True GO-CAM criteria. A True GO-CAM
model is defined as a model where all of the following criteria are met:

- The model is in production status.
- The model has at least two activities that are connected by a causal association, either
  directly or indirectly via shared chemical entities.
- The model has no activities that are disconnected from all other activities in the model.

**Inputs:**

- `--input-dir`: Directory containing GO-CAM model JSON files
- `--output-dir`: Directory to save production True GO-CAM models (required unless using
  `--dry-run`)
- `--pseudo-gocam-output-dir`: Directory to save models that have production status but do not
  otherwise meet True GO-CAM criteria (required unless using `--dry-run`)

**Filtering:**

Models are classified and moved based on these criteria:

- **Status**: Only models with a "production" status are considered.
- **Connectivity**: Models must have at least two activities that are connected by a causal
  association, either
  directly or indirectly via shared chemical entities. Models may not have any activities that are
  completely disconnected from all other activities in the model.

**Outputs:**

- True GO-CAM models are copied to the specified `--output-dir`.
- Models that are production-status but not pathway-like are copied to the
  `--pseudo-gocam-output-dir`.

**Options:**

- `--report-file`: JSON Lines file to write a detailed report of the results
- `--dry-run`: Perform the analysis without copying any files
- `--verbose` / `-v`: Increase verbosity (`-v` for INFO, `-vv` for DEBUG)
- `--limit N`: Limit processing to first N models (0 = no limit)

**Usage:**

```bash
python pipeline/filter_true_gocam_models.py --input-dir /path/to/input --output-dir /path/to/true_models --pseudo-gocam-output-dir /path/to/pseudo_models
```

## add_query_index_to_models.py

Populates the `query_index` field of all models in a directory with pre-computed data for faster
querying and indexing.

**Inputs:**

- `--input-dir`: Directory containing GO-CAM model JSON files
- `--output-dir`: Directory to save indexed GO-CAM model files (required unless using `--dry-run`)

**Options:**

- `--report-file`: JSON Lines file to write a detailed report of the indexing results
- `--dry-run`: Perform indexing without writing output files
- `--go-adapter-descriptor`: OAK adapter descriptor for GO (default: `sqlite:obo:go`)
- `--ncbi-taxon-adapter-descriptor`: OAK adapter descriptor for NCBITaxon (default:
  `sqlite:obo:ncbitaxon`)
- `--goc-groups-yaml`: YAML file defining GOC groups
- `--verbose` / `-v`: Increase verbosity level (`-v` for INFO, `-vv` for DEBUG)
- `--limit N`: Limit processing to first N models (0 = no limit)

**Usage:**

```bash
python pipeline/add_query_index_to_models.py --input-dir /path/to/models --output-dir /path/to/indexed_models
```

## generate_index_files.py

Generates various specialized index files (JSON) mapping attributes like contributors and entities
to model IDs for a directory of GO-CAM models. These files are intended to be used by the GO API.

**Inputs:**

- `--input-dir`: Directory containing indexed GO-CAM model JSON files
- `--output-dir`: Directory to save generated index files (required unless using `--dry-run`)

**Outputs:**

Generates the following JSON index files in the output directory:

- `contributor_index.json`
- `entity_index.json`
- `evidence_index.json`
- `provided_by_index.json`
- `source_index.json`
- `taxon_index.json`

**Options:**

- `--report-file`: JSON Lines file to write a detailed report of the results
- `--dry-run`: Perform processing without writing output files
- `--verbose` / `-v`: Increase verbosity level (`-v` for INFO, `-vv` for DEBUG)
- `--limit N`: Limit processing to first N models (0 = no limit)

**Usage:**

```bash
python pipeline/generate_index_files.py --input-dir /path/to/indexed_models --output-dir /path/to/indices
```

## generate_go_cam_browser_search_docs.py

Generates a single JSON file containing search documents optimized for the GO-CAM Browser.

**Inputs:**

- `--input-dir`: Directory containing indexed GO-CAM model JSON files
- `--output`: File to write the generated search documents to (required unless using `--dry-run`)

**Options:**

- `--report-file`: JSON Lines file to write a detailed report of the results
- `--dry-run`: Perform processing without writing output files
- `--verbose` / `-v`: Increase verbosity level (`-v` for INFO, `-vv` for DEBUG)
- `--limit N`: Limit processing to first N models (0 = no limit)

**Usage:**

```bash
python pipeline/generate_go_cam_browser_search_docs.py --input-dir /path/to/indexed_models --output search.json
```

## generate_log_summary.py

Generates an Excel summary of pipeline run results based on JSONL log files produced by other
pipeline steps.

**Inputs:**

- `--logs-dir`: Directory containing step JSONL report files (e.g., reports generated via
  `--report-file` in other scripts)
- `--output`: File to write the generated Excel summary to (must have `.xlsx` extension)

**Options:**

- `--log-file-extension`: File extension used to find log files (default: `.jsonl`)
- `--verbose` / `-v`: Increase verbosity level (`-v` for INFO, `-vv` for DEBUG)
- `--limit N`: Limit the number of models included in the summary (0 = no limit)

**Usage:**

```bash
python pipeline/generate_log_summary.py --logs-dir /path/to/logs --output summary.xlsx
```

