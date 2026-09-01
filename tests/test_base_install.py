"""Smoke tests for functionality available from the base package installation.

This test is intended to be run in a minimal environment, without any extra dependencies
installed. Typically, tests run with all extra dependencies installed, but this can mask
problems that a user who installs just the base package might encounter.
"""

import importlib
import json

import networkx as nx
import pytest
import yaml
from typer.testing import CliRunner

from gocam import __version__
from gocam.datamodel import Model
from gocam.utils import model_to_digraph
from tests import EXAMPLES_DIR, INPUT_DIR


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def app():
    return importlib.import_module("gocam.cli").app


@pytest.fixture
def release_mock(requests_mock):
    model_id = "663d668500002178"
    model_json = (EXAMPLES_DIR / f"Model-{model_id}.json").read_text()
    requests_mock.get(
        f"https://current.geneontology.org/go-cams/json/{model_id}.json",
        text=model_json,
    )
    return model_id


@pytest.fixture
def local_model():
    with open(EXAMPLES_DIR / "Model-663d668500002178.yaml", "r") as f:
        return Model.model_validate(yaml.safe_load(f))


def test_cli_module_imports():
    assert importlib.import_module("gocam.cli") is not None


def test_version(runner, app):
    result = runner.invoke(app, ["--version"])

    assert result.exit_code == 0
    assert __version__ in result.stdout


def test_fetch_yaml(runner, app, release_mock):
    result = runner.invoke(app, ["fetch", "--format", "yaml", release_mock])

    assert result.exit_code == 0
    parsed_output = yaml.safe_load(result.stdout)
    assert parsed_output["id"] == f"gomodel:{release_mock}"


def test_fetch_json(runner, app, release_mock):
    result = runner.invoke(app, ["fetch", "--format", "json", release_mock])

    assert result.exit_code == 0
    parsed_output = json.loads(result.stdout)
    assert parsed_output["id"] == f"gomodel:{release_mock}"


def test_fetch_from_custom_release(runner, app, requests_mock):
    model_id = "663d668500002178"
    base_url = "https://release.geneontology.org/2026-08-05/go-cams/"
    requests_mock.get(
        f"{base_url}json/{model_id}.json",
        text=(EXAMPLES_DIR / f"Model-{model_id}.json").read_text(),
    )

    result = runner.invoke(
        app,
        ["fetch", "--format", "json", "--base-url", base_url, model_id],
    )

    assert result.exit_code == 0
    assert json.loads(result.stdout)["id"] == f"gomodel:{model_id}"


def test_fetch_without_ids_uses_release_index(runner, app, requests_mock):
    model_id = "663d668500002178"
    base_url = "https://current.geneontology.org/go-cams/"
    requests_mock.get(f"{base_url}index-json/all_index.json", text=f'["{model_id}"]')
    requests_mock.get(
        f"{base_url}json/{model_id}.json",
        text=(EXAMPLES_DIR / f"Model-{model_id}.json").read_text(),
    )

    result = runner.invoke(app, ["fetch", "--format", "json"])

    assert result.exit_code == 0
    assert json.loads(result.stdout)["id"] == f"gomodel:{model_id}"


def test_fetch_as_minerva_warns_and_uses_legacy_endpoint(runner, app, requests_mock):
    model_id = "663d668500002178"
    minerva_object = json.loads((INPUT_DIR / f"minerva-{model_id}.json").read_text())
    requests_mock.get(
        f"https://api.geneontology.org/api/go-cam/{model_id}",
        json=minerva_object,
    )

    result = runner.invoke(app, ["fetch", "--format", "json", "--as-minerva", model_id])

    assert result.exit_code == 0
    assert json.loads(result.stdout)["id"] == f"gomodel:{model_id}"
    assert "deprecated" in result.stderr.lower()


def test_fetch_rejects_base_url_with_as_minerva(runner, app):
    result = runner.invoke(
        app,
        [
            "fetch",
            "--as-minerva",
            "--base-url",
            "https://release.geneontology.org/2026-08-05/go-cams/",
            "663d668500002178",
        ],
    )

    assert result.exit_code == 2
    assert "--base-url" in result.stderr
    assert "--as-minerva" in result.stderr


def test_fetch_all_as_minerva_uses_legacy_index(runner, app, requests_mock):
    model_id = "663d668500002178"
    requests_mock.get(
        "https://go-public.s3.amazonaws.com/files/gocam-models.json",
        json=[{"gocam": f"http://model.geneontology.org/{model_id}"}],
    )
    requests_mock.get(
        f"https://api.geneontology.org/api/go-cam/{model_id}",
        json=json.loads((INPUT_DIR / f"minerva-{model_id}.json").read_text()),
    )

    result = runner.invoke(app, ["fetch", "--format", "json", "--as-minerva"])

    assert result.exit_code == 0
    assert json.loads(result.stdout)["id"] == f"gomodel:{model_id}"
    assert "deprecated" in result.stderr.lower()


def test_load_local_model(local_model):
    assert local_model.id == "gomodel:663d668500002178"


def test_create_networkx_graph(local_model):
    graph = model_to_digraph(local_model)

    assert isinstance(graph, nx.DiGraph)
    assert graph.number_of_nodes() > 0
