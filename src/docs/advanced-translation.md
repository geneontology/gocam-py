# Advanced translation

The base `gocam` package supports fetching and loading GO-CAM models. This page covers
optional translations for CX2 and gene-centered causal networks. Start with
[Getting started](getting-started.md) if you have not yet loaded a model.

## CX2 conversion

Install the CX2 dependencies:

```bash
pip install "gocam[cx2]"
```

The CX2 extra includes `pygraphviz`, which
requires [Graphviz](https://pygraphviz.github.io/documentation/stable/install.html#recommended)
to be installed on the system. Graphviz 2.46 or later is recommended.

Convert a local JSON or YAML model with the CLI:

```bash
gocam convert --output-format cx2 model.yaml > model.cx2.json
```

The input format is inferred from the filename. Run `gocam convert --help` for explicit
input and output formats, file output, graph layout, and optional NDEx upload.

From Python, pass a validated `Model` to `model_to_cx2`:

```python
from gocam.translation.cx2 import model_to_cx2

cx2_document = model_to_cx2(model)
```

Pass `apply_dot_layout=True` to calculate node positions with Graphviz. CLI uploads to
NDEx require the `--ndex-upload` option and the `NDEX_HOST`, `NDEX_USERNAME`, and
`NDEX_PASSWORD` environment variables.

## Gene-to-gene translation

`ModelNetworkTranslator` creates a gene-centered causal network from one or more
validated GO-CAM models. Nodes represent gene products. Edges represent causal
relationships and retain relevant GO terms, evidence references, evidence codes, and
contributors.

Serialize a network in NetworkX node-link JSON format:

```python
from gocam.translation.networkx.model_network_translator import ModelNetworkTranslator

translator = ModelNetworkTranslator()
gene_network_json = translator.translate_models_to_json([model])
```

The serialized document contains:

- `nodes` for gene products, including their identifiers and source model identifiers.
- `edges` for causal relationships, including source and target GO annotations and
  evidence.
- `graph.model_info` for source-model metadata when model information is included.

Translate several models into one combined network by passing them together:

```python
models = [model_a, model_b, model_c]
combined_network_json = translator.translate_models_to_json(models)
```

To work directly with NetworkX rather than serialized JSON:

```python
gene_graph = translator.translate_models([model])
print(gene_graph.number_of_nodes())
print(gene_graph.number_of_edges())
```

## Minerva compatibility

`MinervaWrapper` is available for callers intentionally working with the legacy Minerva
model representation:

```python
from gocam.translation import MinervaWrapper

wrapper = MinervaWrapper()
minerva_object = wrapper.fetch_minerva_object("5b91dbd100002057")
translation = wrapper.translate(minerva_object)
model = translation.result
```

The methods which perform translation internally (`fetch_model()` and `models()`) or use
an outdated index file (`models_ids()`) are deprecated. Use `GoCamClient` for normal
remote retrieval. The CLI's `--as-minerva` option is also deprecated and is retained
temporarily for debugging and migration.
