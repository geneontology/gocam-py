# Getting started

The `gocam` package provides Pydantic models and a command-line interface for working
with Gene Ontology Causal Activity Models (GO-CAMs). See
the [GO-CAM overview](https://geneontology.org/docs/go-annotations/#go-causal-activity-models)
for background on the model and use
the [GO-CAM Browser](https://geneontology.org/go-cam) to explore published models.

## Installation

Install the package from PyPI:

```bash
pip install gocam
```

This installs both the Python package and the `gocam` command.

## Fetch a published model from the command line

The `fetch` command retrieves GO-CAM data from GO Releases. The current release's data
is used by default. YAML is the default output format:

```bash
gocam fetch 5b91dbd100002057 > model.yaml
```

Request JSON with `--format`:

```bash
gocam fetch --format json 5b91dbd100002057 > model.json
```

Use `--base-url` to retrieve models from a dated GO Release:

```bash
gocam fetch \
  --base-url https://release.geneontology.org/2026-08-05/go-cams/ \
  5b91dbd100002057
```

Run `gocam fetch --help` for all fetch options.

## Fetch a published model from Python

The `GoCamClient` class retrieves GO-CAM data from GO Releases. The current release's
data is used by default.

```python
from gocam import GoCamClient

model = GoCamClient().fetch_model("5b91dbd100002057")
print(model.id)
print(model.title)
```

Use `base_url` to retrieve from a dated GO Release:

```python
from gocam import GoCamClient

client = GoCamClient(base_url="https://release.geneontology.org/2026-08-05/go-cams/")
model = client.fetch_model("5b91dbd100002057")
```

## Load a local model

To load a local model, read the file using standard I/O methods and pass the content to
the Pydantic `Model` class for validation.

### JSON

```python
from pathlib import Path

from gocam.datamodel import Model

model = Model.model_validate_json(Path("model.json").read_text())
print(model.id)
```

### YAML

```python
from pathlib import Path

import yaml

from gocam.datamodel import Model

with Path("model.yaml").open() as stream:
    model = Model.model_validate(yaml.safe_load(stream))

print(model.id)
```

## Work with a model

The validated object exposes model fields as Python attributes:

```python
print(model.title)
print(model.taxon)
print(len(model.activities or []))
```

Convert it back to JSON-compatible Python values with Pydantic:

```python
model_data = model.model_dump(mode="json", exclude_none=True)
```

See the [schema reference](index.md) for all model classes and fields. For CX2,
gene-to-gene, and NetworkX workflows, continue
to [Advanced translation](advanced-translation.md).
