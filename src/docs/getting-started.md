# Getting started

The `gocam` package provides Pydantic models and a command-line interface for working with Gene
Ontology Causal Activity Models (GO-CAMs). See the [GO-CAM overview](https://geneontology.org/docs/go-annotations/#go-causal-activity-models) for background
on the model and use the [GO-CAM Browser](https://geneontology.org/go-cam) to explore published models.

## Installation

Install the package from PyPI:

```bash
pip install gocam
```

This installs both the Python package and the `gocam` command. Optional dependencies for CX2 and
gene-network workflows are covered in [Advanced translation](advanced-translation.md).

## Fetch a published model from the command line

The `fetch` command retrieves a published model and translates it into the package's GO-CAM data
model. YAML is the default output format:

```bash
gocam fetch 5b91dbd100002057 > model.yaml
```

Request JSON with `--format`:

```bash
gocam fetch --format json 5b91dbd100002057 > model.json
```

Run `gocam fetch --help` for all fetch options.

## Fetch a published model from Python

`MinervaWrapper` fetches the published Minerva representation and translates it into a validated
`Model`:

```python
from gocam.translation import MinervaWrapper

model = MinervaWrapper().fetch_model("5b91dbd100002057")
print(model.id)
print(model.title)
```

## Load a local model

Loading a local file validates it against the Pydantic model included in the installed `gocam`
version.

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

See the [schema reference](index.md) for all model classes and fields. For CX2, gene-to-gene, and
NetworkX workflows, continue to [Advanced translation](advanced-translation.md).
