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
def api_mock(requests_mock):
    gocam_id = "5b91dbd100002057"
    with open(INPUT_DIR / f"minerva-{gocam_id}.json", "r") as f:
        minerva_object = json.load(f)
    requests_mock.get(
        f"https://api.geneontology.org/api/go-cam/{gocam_id}", json=minerva_object
    )


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


def test_fetch_yaml(runner, app, api_mock):
    result = runner.invoke(app, ["fetch", "--format", "yaml", "5b91dbd100002057"])

    assert result.exit_code == 0
    parsed_output = yaml.safe_load(result.stdout)
    assert parsed_output["id"] == "gomodel:5b91dbd100002057"


def test_fetch_json(runner, app, api_mock):
    result = runner.invoke(app, ["fetch", "--format", "json", "5b91dbd100002057"])

    assert result.exit_code == 0
    parsed_output = json.loads(result.stdout)
    assert parsed_output["id"] == "gomodel:5b91dbd100002057"


def test_load_local_model(local_model):
    assert local_model.id == "gomodel:663d668500002178"


def test_create_networkx_graph(local_model):
    graph = model_to_digraph(local_model)

    assert isinstance(graph, nx.DiGraph)
    assert graph.number_of_nodes() > 0
