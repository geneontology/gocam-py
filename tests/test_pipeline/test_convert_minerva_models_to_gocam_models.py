import json
import shutil
from pathlib import Path

import pytest
from _common import (
    ErrorReason,
    ErrorResult,
    FilteredResult,
    FilterReason,
    SuccessResult,
)
from convert_minerva_models_to_gocam_models import (
    app,
    classify_model,
    process_minerva_model_file,
)
from typer.testing import CliRunner

from gocam.datamodel import (
    Activity,
    CausalAssociation,
    EnabledByGeneProductAssociation,
    EvidenceItem,
    Model,
    ModelStateEnum,
)
from gocam.translation import TranslationWarning, WarningType
from gocam.vocabulary import Relation
from tests import INPUT_DIR

TRUE_GOCAM = "minerva-5b91dbd100002057"
PSEUDO_GOCAM_DISCONNECTED = "minerva-665912ed00002626"
EXCLUDED_NO_ACTIVITY_EDGE = "minerva-YeastPathways_LYSDEGII-PWY"
EXCLUDED_NOT_PRODUCTION = "minerva-nested-localization"


def activity(
    activity_id: str, *, downstream: str | None = None, evidence: bool = True
) -> Activity:
    return Activity(
        id=activity_id,
        enabled_by=EnabledByGeneProductAssociation(
            term=f"UniProtKB:{activity_id}",
            evidence=(
                [EvidenceItem(term="ECO:0000314", reference="PMID:1")]
                if evidence
                else None
            ),
        ),
        causal_associations=(
            [
                CausalAssociation(
                    predicate=Relation.CAUSALLY_UPSTREAM_OF,
                    downstream_activity=downstream,
                )
            ]
            if downstream
            else None
        ),
    )


def make_model(
    *activities: Activity, status: ModelStateEnum | None = ModelStateEnum.production
) -> Model:
    return Model(
        id="gomodel:test",
        title="Test model",
        status=status,
        activities=list(activities),
    )


def connected(**kwargs) -> Model:
    return make_model(activity("a", downstream="b"), activity("b"), **kwargs)


def disconnected() -> Model:
    return make_model(activity("a", downstream="b"), activity("b"), activity("c"))


def no_evidence() -> Model:
    return make_model(
        activity("a", downstream="b", evidence=False), activity("b", evidence=False)
    )


COMPLEMENT = TranslationWarning(
    type=WarningType.INDIVIDUAL_HAS_COMPLEMENT_TYPE, message="complement"
)
MOLECULAR_EVENT = TranslationWarning(
    type=WarningType.MOLECULAR_EVENT_INDIVIDUAL, message="molecular event"
)
UNRELATED = TranslationWarning(type=WarningType.MISSING_TERM, message="missing term")


@pytest.mark.parametrize(
    ("model", "warnings", "expected_reason"),
    [
        pytest.param(connected(), [], None, id="true-gocam"),
        pytest.param(connected(), [UNRELATED], None, id="unrelated-warning"),
        pytest.param(
            connected(status=ModelStateEnum.development),
            [COMPLEMENT],
            FilterReason.NOT_PRODUCTION_MODEL,
            id="not-production-beats-complement",
        ),
        pytest.param(
            connected(status=None),
            [],
            FilterReason.NOT_PRODUCTION_MODEL,
            id="no-status",
        ),
        pytest.param(
            make_model(activity("a")),
            [MOLECULAR_EVENT],
            FilterReason.NO_ACTIVITY_EDGE,
            id="no-edge-beats-molecular-event",
        ),
        pytest.param(
            make_model(),
            [],
            FilterReason.NO_ACTIVITY_EDGE,
            id="no-activities",
        ),
        pytest.param(
            connected(),
            [MOLECULAR_EVENT, COMPLEMENT],
            FilterReason.USES_COMPLEMENT,
            id="complement-beats-molecular-event",
        ),
        pytest.param(
            disconnected(),
            [MOLECULAR_EVENT],
            FilterReason.USES_MOLECULAR_EVENT,
            id="molecular-event-beats-disconnected",
        ),
        pytest.param(
            disconnected(),
            [],
            FilterReason.DISCONNECTED_ACTIVITY,
            id="disconnected",
        ),
        pytest.param(no_evidence(), [], FilterReason.NO_EVIDENCE, id="no-evidence"),
    ],
)
def test_classify_model(
    model: Model,
    warnings: list[TranslationWarning],
    expected_reason: FilterReason | None,
):
    result = classify_model(model, warnings)

    if expected_reason is None:
        assert result == SuccessResult()
    else:
        assert result == FilteredResult(reason=expected_reason)


def copy_fixtures(destination: Path, *stems: str) -> None:
    destination.mkdir()
    for stem in stems:
        shutil.copy(INPUT_DIR / f"{stem}.json", destination)


def read_log(log_file: Path) -> tuple[dict, dict[str, dict]]:
    header, *entries = (json.loads(line) for line in log_file.read_text().splitlines())
    return header, {entry["model_id"]: entry for entry in entries}


def json_names(directory: Path) -> set[str]:
    return {path.name for path in directory.glob("*.json")}


def test_cli_sorts_models_into_output_dirs(tmp_path: Path):
    input_dir = tmp_path / "input"
    copy_fixtures(
        input_dir,
        TRUE_GOCAM,
        PSEUDO_GOCAM_DISCONNECTED,
        EXCLUDED_NO_ACTIVITY_EDGE,
        EXCLUDED_NOT_PRODUCTION,
    )
    (input_dir / "broken.json").write_text("{")
    true_dir = tmp_path / "true"
    pseudo_dir = tmp_path / "pseudo"
    true_dir.mkdir()
    pseudo_dir.mkdir()
    log_file = tmp_path / "convert.jsonl"

    result = CliRunner().invoke(
        app,
        [
            "--input-dir",
            str(input_dir),
            "--output-dir",
            str(true_dir),
            "--pseudo-gocam-output-dir",
            str(pseudo_dir),
            "--log-file",
            str(log_file),
        ],
    )

    assert result.exit_code == 0, result.output
    assert json_names(true_dir) == {f"{TRUE_GOCAM}.json"}
    assert json_names(pseudo_dir) == {f"{PSEUDO_GOCAM_DISCONNECTED}.json"}
    Model.model_validate_json(
        (pseudo_dir / f"{PSEUDO_GOCAM_DISCONNECTED}.json").read_text()
    )

    header, entries = read_log(log_file)
    assert header == {
        "record_type": "pipeline_step_log",
        "format_version": 2,
        "step": "convert",
    }
    assert entries[TRUE_GOCAM]["status"] == "success"
    assert "category" not in entries[TRUE_GOCAM]
    assert {
        model_id: (entry["status"], entry.get("reason"), entry.get("category"))
        for model_id, entry in entries.items()
        if model_id != TRUE_GOCAM
    } == {
        PSEUDO_GOCAM_DISCONNECTED: (
            "filtered",
            FilterReason.DISCONNECTED_ACTIVITY.value,
            "pseudo-gocam",
        ),
        EXCLUDED_NO_ACTIVITY_EDGE: (
            "filtered",
            FilterReason.NO_ACTIVITY_EDGE.value,
            "excluded",
        ),
        EXCLUDED_NOT_PRODUCTION: (
            "filtered",
            FilterReason.NOT_PRODUCTION_MODEL.value,
            "excluded",
        ),
        "broken": ("error", ErrorReason.READ_ERROR.value, None),
    }
    assert entries[PSEUDO_GOCAM_DISCONNECTED]["meta"]["status"] == "production"


def test_cli_dry_run_writes_log_only(tmp_path: Path):
    input_dir = tmp_path / "input"
    copy_fixtures(input_dir, TRUE_GOCAM, PSEUDO_GOCAM_DISCONNECTED)
    log_file = tmp_path / "convert.jsonl"

    result = CliRunner().invoke(
        app,
        ["--input-dir", str(input_dir), "--log-file", str(log_file), "--dry-run"],
    )

    assert result.exit_code == 0, result.output
    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "convert.jsonl",
        "input",
    ]
    _, entries = read_log(log_file)
    assert entries[TRUE_GOCAM]["status"] == "success"
    assert entries[PSEUDO_GOCAM_DISCONNECTED]["category"] == "pseudo-gocam"


def test_cli_requires_pseudo_gocam_output_dir(tmp_path: Path):
    input_dir = tmp_path / "input"
    copy_fixtures(input_dir, TRUE_GOCAM)
    true_dir = tmp_path / "true"
    true_dir.mkdir()

    result = CliRunner().invoke(
        app, ["--input-dir", str(input_dir), "--output-dir", str(true_dir)]
    )

    assert result.exit_code == 2
    assert json_names(true_dir) == set()


def test_pseudo_gocam_write_error_keeps_meta(tmp_path: Path):
    true_dir = tmp_path / "true"
    true_dir.mkdir()
    missing_pseudo_dir = tmp_path / "does-not-exist"

    result = process_minerva_model_file(
        INPUT_DIR / f"{PSEUDO_GOCAM_DISCONNECTED}.json",
        output_dir=true_dir,
        pseudo_gocam_output_dir=missing_pseudo_dir,
    )

    assert isinstance(result, ErrorResult)
    assert result.reason == ErrorReason.WRITE_ERROR
    assert result.meta is not None
    assert result.meta["status"] == "production"


def test_filtered_model_keeps_translation_warnings(tmp_path: Path):
    minerva_model = json.loads(
        (INPUT_DIR / "minerva-5ce58dde00001215.json").read_text()
    )
    for annotation in minerva_model["annotations"]:
        if annotation["key"] == "state":
            annotation["value"] = "development"
    json_file = tmp_path / "development-model.json"
    json_file.write_text(json.dumps(minerva_model))

    result = process_minerva_model_file(
        json_file, output_dir=None, pseudo_gocam_output_dir=None
    )

    assert isinstance(result, FilteredResult)
    assert result.reason == FilterReason.NOT_PRODUCTION_MODEL
    assert len(result.warnings) == 7
    assert result.get_log_entry("development-model")["warnings"] == result.warnings
