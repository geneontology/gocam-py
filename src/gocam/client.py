from collections.abc import Iterator
from dataclasses import dataclass, field
from urllib.parse import quote

import requests
from pydantic import TypeAdapter

from gocam.datamodel import Model

DEFAULT_GO_CAM_BASE_URL = "https://current.geneontology.org/go-cams/"

_MODEL_ID_PREFIXES = ("gomodel:", "gocam:")
_MODEL_IDS_ADAPTER = TypeAdapter(list[str])


@dataclass
class GoCamClient:
    """Client for JSON-serialized standard format models in GO Release data."""

    base_url: str = DEFAULT_GO_CAM_BASE_URL
    """Base URL for GO-CAM model data. Defaults to the current release."""

    session: requests.Session = field(default_factory=requests.Session)
    """Requests session for HTTP requests. Defaults to a new session."""

    timeout: float = 30.0
    """Timeout in seconds for HTTP requests. Defaults to 30 seconds."""

    def fetch_model(self, model_id: str) -> Model:
        """Fetch one GO-CAM model."""

        local_id = self._normalize_model_id(model_id)
        encoded_id = quote(local_id, safe="")
        response = self.session.get(
            self._url(f"json/{encoded_id}.json"), timeout=self.timeout
        )
        response.raise_for_status()
        return Model.model_validate_json(response.text)

    def fetch_model_ids(self) -> list[str]:
        """Fetch the index of all GO-CAM model IDs."""

        response = self.session.get(
            self._url("index-json/all_index.json"), timeout=self.timeout
        )
        response.raise_for_status()
        return _MODEL_IDS_ADAPTER.validate_json(response.text)

    def iter_models(self) -> Iterator[Model]:
        """Lazily fetch models in the order provided by the all-model index."""

        for model_id in self.fetch_model_ids():
            yield self.fetch_model(model_id)

    def _url(self, path: str) -> str:
        return f"{self.base_url.rstrip('/')}/{path}"

    @staticmethod
    def _normalize_model_id(model_id: str) -> str:
        if not model_id:
            raise ValueError(f"Missing GO-CAM ID: {model_id}")
        for prefix in _MODEL_ID_PREFIXES:
            if model_id.startswith(prefix):
                model_id = model_id.removeprefix(prefix)
                break
        if not model_id:
            raise ValueError(f"Missing GO-CAM ID: {model_id}")
        return model_id
