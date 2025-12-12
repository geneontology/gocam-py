# Pipeline Scripts

This directory contains scripts for processing steps in the GOC pipeline.

## convert_minerva_models_to_gocam_models.py

Converts Minerva model JSON files to GO-CAM model format, doing basic filtering to eliminate obvious non-GO-CAM models. The output files are intended for internal QC analysis and further filtering.

**Inputs:**
- `--input-dir`: Directory containing Minerva model JSON files
- `--output-dir`: Directory to save converted GO-CAM models (required unless using `--dry-run`)

**Filtering:**

Models are filtered out and not written if they:
- Have no activity edges (no causal associations between any activities)
- Use complement types

**Outputs:**
- GO-CAM model JSON files (one per input file) in the specified output directory

**Options:**
- `--dry-run`: Perform conversion without writing output files
- `--verbose` / `-v`: Increase verbosity (use multiple times for more detail: `-v` for INFO, `-vv` for DEBUG)
- `--limit N`: Limit processing to first N models (0 = no limit)

**Usage:**
```bash
python convert_minerva_models_to_gocam_models.py --input-dir /path/to/minerva/models --output-dir /path/to/gocam/models
```

