# Pipeline Scripts

This directory contains scripts for processing steps in the GOC pipeline.

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

Filters a collection of models based on whether they meet the "True GO-CAM" criteria. A True GO-CAM model is defined as a production model that is pathway-like (containing at least three activities connected in sequence).

**Inputs:**

- `--input-dir`: Directory containing GO-CAM model JSON files
- `--output-dir`: Directory to save production True GO-CAM models (required unless using `--dry-run`)
- `--pseudo-gocam-output-dir`: Directory to save production pseudo-GO-CAM models (required unless using `--dry-run`)

**Filtering:**

Models are classified and moved based on these criteria:

- **Status**: Only models with a "production" status are considered.
- **Connectivity**: Models must be "pathway-like", meaning they have a path of at least three activities connected via causal associations or shared chemical inputs/outputs.

**Outputs:**

- True GO-CAM models are copied to the specified `--output-dir`.
- Models that are production-status but not pathway-like are copied to the `--pseudo-gocam-output-dir`.

**Options:**

- `--report-file`: JSON Lines file to write a detailed report of the results
- `--dry-run`: Perform the analysis without copying any files
- `--verbose` / `-v`: Increase verbosity (`-v` for INFO, `-vv` for DEBUG)
- `--limit N`: Limit processing to first N models (0 = no limit)

**Usage:**

```bash
python pipeline/filter_true_gocam_models.py --input-dir /path/to/input --output-dir /path/to/true_models --pseudo-gocam-output-dir /path/to/pseudo_models
```

