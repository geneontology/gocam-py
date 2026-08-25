"""Tests for the command-line interface (CLI) of the gocam package.

Note: important basic tests of the CLI are in tests/test_base_install.py so that they
can be run separately in a minimal environment in CI.
"""

import json

import pytest
from typer.testing import CliRunner

from gocam.cli import app
from tests import EXAMPLES_DIR


@pytest.fixture
def runner():
    return CliRunner()


@pytest.mark.parametrize("format", ["json", "yaml"])
def test_convert_to_cx2_from_file(runner, format):
    result = runner.invoke(
        app,
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
            app, ["convert", "-O", "cx2", "-I", format], input=f.read()
        )
    assert result.exit_code == 0
    cx2 = json.loads(result.stdout)
    assert isinstance(cx2, list)


def test_convert_to_cx2_to_file(runner, tmp_path):
    output_path = tmp_path / "test.cx2"
    result = runner.invoke(
        app,
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
