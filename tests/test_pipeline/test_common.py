from _common import (
    PIPELINE_LOG_FORMAT_VERSION,
    FilterCategory,
    FilteredResult,
    FilterReason,
)

EXPECTED_CATEGORIES = {
    FilterReason.NOT_PRODUCTION_MODEL: FilterCategory.EXCLUDED,
    FilterReason.NO_ACTIVITY_EDGE: FilterCategory.EXCLUDED,
    FilterReason.DISCONNECTED_ACTIVITY: FilterCategory.PSEUDO_GOCAM,
    FilterReason.NO_EVIDENCE: FilterCategory.PSEUDO_GOCAM,
    FilterReason.USES_COMPLEMENT: FilterCategory.PSEUDO_GOCAM,
    FilterReason.USES_MOLECULAR_EVENT: FilterCategory.PSEUDO_GOCAM,
}


def test_every_filter_reason_has_expected_category():
    assert {reason: reason.category for reason in FilterReason} == EXPECTED_CATEGORIES


def test_filtered_result_log_entry_includes_category():
    result = FilteredResult(reason=FilterReason.USES_MOLECULAR_EVENT)

    assert result.get_log_entry("gomodel:123") == {
        "model_id": "gomodel:123",
        "status": "filtered",
        "reason": "Model uses molecular event",
        "category": "pseudo-gocam",
    }


def test_pipeline_log_format_version():
    assert PIPELINE_LOG_FORMAT_VERSION == 2


def test_filtered_result_log_entry_includes_warnings():
    warning = {"type": "MISSING_TERM", "message": "missing term", "entity_id": None}
    result = FilteredResult(reason=FilterReason.NO_EVIDENCE, warnings=[warning])

    assert result.get_log_entry("gomodel:123")["warnings"] == [warning]
