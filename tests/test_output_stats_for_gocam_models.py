"""Tests for ``pipeline/output_stats_for_gocam_models.py``. The stats script
lives under ``pipeline/`` (excluded from the package), so it is loaded by file
path via importlib. Add tests for any update to that module here."""

import importlib.util
from pathlib import Path

from typer.testing import CliRunner

from gocam.datamodel import (
    Activity,
    CausalAssociation,
    EnabledByProteinComplexAssociation,
    Model,
    ModelStateEnum,
    MolecularFunctionAssociation,
)

PIPELINE_SCRIPT = (
    Path(__file__).parent.parent / "pipeline" / "output_stats_for_gocam_models.py"
)


def _load_stats_module():
    spec = importlib.util.spec_from_file_location(
        "output_stats_for_gocam_models", PIPELINE_SCRIPT
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


stats = _load_stats_module()


def _sample_model() -> Model:
    """A production model with a protein-complex activity and a
    constitutively-upstream relation, so both record types are non-empty."""
    return Model(
        id="gomodel:TEST",
        title="Test model",
        status=ModelStateEnum.production,
        activities=[
            Activity(
                id="a1",
                enabled_by=EnabledByProteinComplexAssociation(term="GO:0000151"),
                molecular_function=MolecularFunctionAssociation(term="GO:0003674"),
                causal_associations=[
                    CausalAssociation(predicate="RO:0012009", downstream_activity="a2")
                ],
            ),
            Activity(
                id="a2",
                enabled_by=EnabledByProteinComplexAssociation(term="GO:0000152"),
            ),
        ],
    )


def _write_model_json(dir_path: Path, model: Model) -> None:
    dir_path.mkdir(parents=True, exist_ok=True)
    (dir_path / "model.json").write_text(model.model_dump_json())


def test_write_tsv_aggregate_files_skips_empty(tmp_path):
    out = tmp_path / "all_models"
    stats.write_tsv_aggregate_files(out, [], [])
    assert not out.exists()  # no directory, no files


def test_write_tsv_aggregate_files_writes_only_nonempty(tmp_path):
    out = tmp_path / "all_models"
    pc = [stats.ProteinComplexActivityInfo(model_id="m", activity_id="a")]
    stats.write_tsv_aggregate_files(out, pc, [])
    assert (out / "aggregate_protein_complex.json").exists()
    assert not (out / "aggregate_constitutively_upstream.json").exists()


def test_all_models_requires_tsv_dir(tmp_path):
    runner = CliRunner()
    result = runner.invoke(
        stats.app,
        ["--input-dir", str(tmp_path), "--dry-run", "--all-models"],
    )
    assert result.exit_code != 0
    assert "all-models-tsv-dir" in result.output


def test_default_no_all_models_dir_even_with_records(tmp_path):
    input_dir = tmp_path / "in"
    _write_model_json(input_dir, _sample_model())
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    all_models_dir = tmp_path / "all_models"
    runner = CliRunner()
    # --all-models NOT passed => default off; the dir must be ignored.
    result = runner.invoke(
        stats.app,
        [
            "--input-dir",
            str(input_dir),
            "--output-dir",
            str(output_dir),
            "--all-models-tsv-dir",
            str(all_models_dir),
        ],
    )
    assert result.exit_code == 0, result.output
    assert not all_models_dir.exists()
    # The production protein-complex / constitutively-upstream files are gated
    # on --all-models too, so they must not appear in the main output dir.
    assert not (output_dir / "aggregate_protein_complex.json").exists()
    assert not (output_dir / "aggregate_constitutively_upstream.json").exists()


def test_all_models_on_writes_record_files(tmp_path):
    input_dir = tmp_path / "in"
    _write_model_json(input_dir, _sample_model())
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    all_models_dir = tmp_path / "all_models"
    runner = CliRunner()
    result = runner.invoke(
        stats.app,
        [
            "--input-dir",
            str(input_dir),
            "--output-dir",
            str(output_dir),
            "--all-models",
            "--all-models-tsv-dir",
            str(all_models_dir),
        ],
    )
    assert result.exit_code == 0, result.output
    assert (all_models_dir / "aggregate_protein_complex.json").exists()
    assert (all_models_dir / "aggregate_constitutively_upstream.json").exists()
    # With --all-models on, the production (status-filtered) versions are also
    # written into the main output dir.
    assert (output_dir / "aggregate_protein_complex.json").exists()
    assert (output_dir / "aggregate_constitutively_upstream.json").exists()
