# gocam

`gocam` is a Python data model and command-line interface for working with Gene Ontology Causal
Activity Models (GO-CAMs). It provides Pydantic models for validated GO-CAM data, tools for fetching
published models, and optional translations to formats such as CX2 and gene-centered networks.

Learn more about
[GO-CAMs](https://geneontology.org/docs/go-annotations/#go-causal-activity-models) or explore
published models in the [GO-CAM browser](https://geneontology.org/go-cam).

## Installation

Install `gocam` from PyPI:

```bash
pip install gocam
```

This installs both the Python package and the `gocam` command.

## Quick start

Fetch a published model as YAML from the command line:

```bash
gocam fetch 5b91dbd100002057 > model.yaml
```

Fetch and translate a published model from Python:

```python
from gocam.translation import MinervaWrapper

model = MinervaWrapper().fetch_model("5b91dbd100002057")
print(model.title)
```

Load and validate a local JSON model:

```python
from pathlib import Path

from gocam.datamodel import Model

model = Model.model_validate_json(Path("model.json").read_text())
print(model.id)
```

See the [getting-started guide](https://geneontology.github.io/gocam-py/getting-started/) for JSON
output, local YAML loading, serialization, and more complete examples.

## Documentation

- [Getting started](https://geneontology.github.io/gocam-py/getting-started/)
- [Generated schema reference](https://geneontology.github.io/gocam-py/)
- [Advanced translation](https://geneontology.github.io/gocam-py/advanced-translation/)
- [Versioned schema artifacts](https://github.com/geneontology/gocam-py/releases)

## Development

The source schema is `src/gocam/schema/gocam.yaml`. This repository also contains reusable Python
code and internal scripts used by the GO Release pipeline. One-off and exploratory analyses belong
in the [`gocam-analysis` repository](https://github.com/geneontology/gocam-analysis).

Install [uv](https://docs.astral.sh/uv/), clone the repository, and create the complete development
environment:

```bash
git clone https://github.com/geneontology/gocam-py.git
cd gocam-py
uv sync --all-extras
```

CX2-related development requires
[Graphviz](https://pygraphviz.github.io/documentation/stable/install.html#recommended) to be
installed on the system before syncing the optional dependencies.

Run the standard project checks:

```bash
make test
make lint-python
make type-check
```

Regenerate schema-derived outputs when intentionally changing the schema:

```bash
make src/gocam/datamodel/gocam.py
make gen-project
make gendoc
```

The generated `project/` directory is disposable local output. Versioned generated schema artifacts
are distributed through [GitHub Releases](https://github.com/geneontology/gocam-py/releases).
