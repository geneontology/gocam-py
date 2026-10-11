from pathlib import Path
from unittest.mock import Mock

import pytest
import requests
from pydantic import ValidationError

from gocam import GoCamClient
from tests import EXAMPLES_DIR


@pytest.fixture
def model_json() -> str:
    return Path(EXAMPLES_DIR / "Model-663d668500002178.json").read_text()


@pytest.mark.parametrize(
    "model_id",
    [
        "663d668500002178",
        "gomodel:663d668500002178",
        "gocam:663d668500002178",
    ],
)
def test_fetch_model_normalizes_id_and_validates_model(
    requests_mock, model_json, model_id
):
    requests_mock.get(
        "https://current.geneontology.org/go-cams/json/663d668500002178.json",
        text=model_json,
    )

    model = GoCamClient().fetch_model(model_id)

    assert model.id == "gomodel:663d668500002178"


def test_fetch_model_uses_custom_release_base_url(requests_mock, model_json):
    requests_mock.get(
        "https://release.geneontology.org/2026-08-05/go-cams/json/"
        "663d668500002178.json",
        text=model_json,
    )
    client = GoCamClient(base_url="https://release.geneontology.org/2026-08-05/go-cams")

    model = client.fetch_model("663d668500002178")

    assert model.id == "gomodel:663d668500002178"


@pytest.mark.parametrize("model_id", ["", "gomodel:", "gocam:"])
def test_fetch_model_rejects_empty_id(model_id):
    with pytest.raises(ValueError, match="Missing GO-CAM ID"):
        GoCamClient().fetch_model(model_id)


def test_fetch_model_encodes_id_as_one_path_segment(model_json):
    response = Mock(spec=requests.Response)
    response.text = model_json
    session = Mock(spec=requests.Session)
    session.get.return_value = response

    GoCamClient(session=session).fetch_model("model/child")

    session.get.assert_called_once_with(
        "https://current.geneontology.org/go-cams/json/model%2Fchild.json",
        timeout=30.0,
    )


def test_fetch_model_propagates_http_error(requests_mock):
    requests_mock.get(
        "https://current.geneontology.org/go-cams/json/missing.json",
        status_code=404,
    )

    with pytest.raises(requests.HTTPError):
        GoCamClient().fetch_model("missing")


def test_fetch_model_propagates_transport_error():
    session = Mock(spec=requests.Session)
    session.get.side_effect = requests.ConnectionError("connection failed")

    with pytest.raises(requests.ConnectionError, match="connection failed"):
        GoCamClient(session=session).fetch_model("unreachable")


@pytest.mark.parametrize("response_text", ["not-json", "[]"])
def test_fetch_model_rejects_invalid_standard_form(requests_mock, response_text):
    requests_mock.get(
        "https://current.geneontology.org/go-cams/json/invalid.json",
        text=response_text,
    )

    with pytest.raises(ValidationError):
        GoCamClient().fetch_model("invalid")


def test_fetch_model_forwards_timeout(model_json):
    response = Mock(spec=requests.Response)
    response.text = model_json
    session = Mock(spec=requests.Session)
    session.get.return_value = response
    client = GoCamClient(session=session, timeout=12.5)

    client.fetch_model("663d668500002178")

    session.get.assert_called_once_with(
        "https://current.geneontology.org/go-cams/json/663d668500002178.json",
        timeout=12.5,
    )
    response.raise_for_status.assert_called_once_with()


def test_fetch_model_ids_validates_string_list(requests_mock):
    requests_mock.get(
        "https://current.geneontology.org/go-cams/index-json/all_index.json",
        text='["model-b", "model-a"]',
    )

    assert GoCamClient().fetch_model_ids() == ["model-b", "model-a"]


@pytest.mark.parametrize("index_json", ["not-json", '{"ALL": []}', '["model-a", 7]'])
def test_fetch_model_ids_rejects_invalid_index(requests_mock, index_json):
    requests_mock.get(
        "https://current.geneontology.org/go-cams/index-json/all_index.json",
        text=index_json,
    )

    with pytest.raises(ValidationError):
        GoCamClient().fetch_model_ids()


def test_iter_models_is_lazy_and_preserves_index_order(requests_mock, model_json):
    index_url = "https://current.geneontology.org/go-cams/index-json/all_index.json"
    first_url = "https://current.geneontology.org/go-cams/json/model-a.json"
    second_url = "https://current.geneontology.org/go-cams/json/model-b.json"
    requests_mock.get(index_url, text='["model-a", "model-b"]')
    requests_mock.get(first_url, text=model_json)
    requests_mock.get(second_url, text=model_json)

    models = GoCamClient().iter_models()
    assert requests_mock.request_history == []

    next(models)
    assert [request.url for request in requests_mock.request_history] == [
        index_url,
        first_url,
    ]

    next(models)
    assert [request.url for request in requests_mock.request_history] == [
        index_url,
        first_url,
        second_url,
    ]
