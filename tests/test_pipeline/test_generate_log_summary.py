import json
from pathlib import Path

import pytest
from _common import (
    PIPELINE_LOG_FORMAT_VERSION,
    PIPELINE_LOG_RECORD_TYPE,
    ErrorReason,
    FilterReason,
    PipelineStep,
)
from generate_log_summary import (
    ModelSummary,
    build_group_reports,
    discover_step_logs,
    summarize_model_entries,
)


def model_summary(
    model_id: str,
    *,
    status: str | None = "production",
    convert_status: str | None = None,
    convert_category: str | None = None,
    convert_reason: str | None = None,
) -> ModelSummary:
    return ModelSummary(
        model_id=model_id,
        model_status=status,
        convert_status=convert_status,
        convert_category=convert_category,
        convert_reason=convert_reason,
    )


def test_build_group_reports_selects_models_by_convert_outcome():
    summaries = [
        model_summary("true-gocam", convert_status="success"),
        model_summary(
            "pseudo-gocam",
            convert_status="filtered",
            convert_category="pseudo-gocam",
            convert_reason=FilterReason.USES_MOLECULAR_EVENT.value,
        ),
        model_summary(
            "write-error",
            convert_status="error",
            convert_reason=ErrorReason.WRITE_ERROR.value,
        ),
        model_summary(
            "excluded",
            convert_status="filtered",
            convert_category="excluded",
            convert_reason=FilterReason.NO_ACTIVITY_EDGE.value,
        ),
        model_summary(
            "conversion-error",
            convert_status="error",
            convert_reason=ErrorReason.CONVERSION_ERROR.value,
        ),
        model_summary(
            "read-error",
            status=None,
            convert_status="error",
            convert_reason=ErrorReason.READ_ERROR.value,
        ),
        model_summary(
            "not-production",
            status="development",
            convert_status="filtered",
            convert_category="excluded",
            convert_reason=FilterReason.NOT_PRODUCTION_MODEL.value,
        ),
        model_summary("no-convert-entry"),
    ]

    reports = build_group_reports(summaries)

    assert {model.model_id for report in reports for model in report.models} == {
        "true-gocam",
        "pseudo-gocam",
        "write-error",
    }


def test_summarize_model_entries_records_convert_outcome():
    summary = summarize_model_entries(
        "m1",
        [
            (
                PipelineStep.CONVERT,
                {
                    "model_id": "m1",
                    "status": "filtered",
                    "reason": FilterReason.USES_COMPLEMENT.value,
                    "category": "pseudo-gocam",
                    "meta": {"status": "production", "title": "Model 1"},
                },
            ),
        ],
    )

    assert summary.convert_status == "filtered"
    assert summary.convert_category == "pseudo-gocam"
    assert summary.convert_reason == FilterReason.USES_COMPLEMENT.value
    assert summary.model_status == "production"
    assert summary.pipeline_status == "filtered"
    assert summary.pipeline_status_details == FilterReason.USES_COMPLEMENT.value


def test_summarize_model_entries_ignores_other_steps_for_convert_outcome():
    summary = summarize_model_entries(
        "m1",
        [
            (
                PipelineStep.QUERY_INDEX,
                {
                    "model_id": "m1",
                    "status": "success",
                    "meta": {"status": "production"},
                },
            ),
        ],
    )

    assert summary.convert_status is None
    assert summary.convert_category is None
    assert summary.convert_reason is None
    assert build_group_reports([summary]) == []


def write_log_header(path: Path, *, format_version: int, step: str) -> None:
    header = {
        "record_type": PIPELINE_LOG_RECORD_TYPE,
        "format_version": format_version,
        "step": step,
    }
    path.write_text(json.dumps(header) + "\n")


def test_discover_step_logs_rejects_format_version_1(tmp_path: Path):
    write_log_header(tmp_path / "01-convert.jsonl", format_version=1, step="convert")

    with pytest.raises(ValueError, match="unsupported pipeline log format version 1"):
        discover_step_logs(tmp_path)


def test_discover_step_logs_rejects_filter_step(tmp_path: Path):
    write_log_header(
        tmp_path / "02-filter.jsonl",
        format_version=PIPELINE_LOG_FORMAT_VERSION,
        step="filter",
    )

    with pytest.raises(ValueError, match="unknown pipeline step 'filter'"):
        discover_step_logs(tmp_path)
