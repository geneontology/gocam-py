import json

import pytest
import yaml
from click.testing import CliRunner

from gocam import __version__
from gocam.cli import fetch, cli
from tests import INPUT_DIR


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def api_mock(requests_mock):
    gocam_id = "5b91dbd100002057"
    with open(INPUT_DIR / f"minerva-{gocam_id}.json", "r") as f:
        minerva_object = json.load(f)
    requests_mock.get(
        f"https://api.geneontology.org/api/go-cam/{gocam_id}", json=minerva_object
    )


def test_fetch_yaml(runner, api_mock):
    result = runner.invoke(fetch, ["--format", "yaml", "5b91dbd100002057"])
    assert result.exit_code == 0

    parsed_output = yaml.safe_load(result.output)
    assert parsed_output["id"] == "gomodel:5b91dbd100002057"


def test_fetch_json(runner, api_mock):
    result = runner.invoke(fetch, ["--format", "json", "5b91dbd100002057"])
    assert result.exit_code == 0

    parsed_output = json.loads(result.output)
    assert parsed_output["id"] == "gomodel:5b91dbd100002057"


def test_version(runner):
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output
