import json

import pytest
import yaml
from click.testing import CliRunner

from gocam import __version__
from gocam.cli import cli
from tests import EXAMPLES_DIR, INPUT_DIR


@pytest.fixture
def runner():
    return CliRunner(mix_stderr=False)


@pytest.fixture
def api_mock(requests_mock):
    gocam_id = "5b91dbd100002057"
    with open(INPUT_DIR / f"minerva-{gocam_id}.json", "r") as f:
        minerva_object = json.load(f)
    requests_mock.get(
        f"https://api.geneontology.org/api/go-cam/{gocam_id}", json=minerva_object
    )


def test_fetch_yaml(runner, api_mock):
    result = runner.invoke(cli, ["fetch", "--format", "yaml", "5b91dbd100002057"])
    assert result.exit_code == 0

    parsed_output = yaml.safe_load(result.stdout)
    assert parsed_output["id"] == "gomodel:5b91dbd100002057"


def test_fetch_json(runner, api_mock):
    result = runner.invoke(cli, ["fetch", "--format", "json", "5b91dbd100002057"])
    assert result.exit_code == 0

    parsed_output = json.loads(result.stdout)
    assert parsed_output["id"] == "gomodel:5b91dbd100002057"


def test_version(runner):
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.stdout


@pytest.mark.parametrize("format", ["json", "yaml"])
def test_convert_to_cx2_from_file(runner, format):
    result = runner.invoke(
        cli,
        [
            "convert",
            "-O",
            "cx2",
            str(EXAMPLES_DIR / f"Model-663d668500002178.{format}"),
        ],
    )
    assert result.exit_code == 0
    cx2 = json.loads(result.stdout)
    assert isinstance(cx2, list)


@pytest.mark.parametrize("format", ["json", "yaml"])
def test_convert_to_cx2_from_stdin(runner, format):
    with open(EXAMPLES_DIR / f"Model-663d668500002178.{format}") as f:
        result = runner.invoke(
            cli, ["convert", "-O", "cx2", "-I", format], input=f.read()
        )
    assert result.exit_code == 0
    cx2 = json.loads(result.stdout)
    assert isinstance(cx2, list)


def test_convert_to_cx2_to_file(runner, tmp_path):
    output_path = tmp_path / "test.cx2"
    result = runner.invoke(
        cli,
        [
            "convert",
            "-O",
            "cx2",
            "-o",
            str(output_path),
            str(EXAMPLES_DIR / "Model-663d668500002178.yaml"),
        ],
    )

    assert result.exit_code == 0
    assert output_path.exists()
    with open(output_path) as f:
        cx2 = json.load(f)
    assert isinstance(cx2, list)
