from dataclasses import dataclass
from typing import Union, List, Dict, Optional, Any

import logging
from gocam.datamodel import Model

logger = logging.getLogger(__name__)

FIELD_VALUE = Union[str, int, float, bool, None, List[str]]
ROW = Dict[str, FIELD_VALUE]


@dataclass
class Flattener:
    """
    A class to flatten a GO-CAM model into a single row representation.

    Each field in the row is either:

    - an atomic value (str, int, float, bool, None)
    - a list of strings (for fields that are lists in the model, e.g. terms)
    """
    fields: Optional[List[str]] = None

    def flatten(self, model: Model) -> ROW:
        """
        Flatten a GO-CAM model into a single row representation.

        Args:
            model (Model): The GO-CAM model to flatten.

        Returns:
            ROW: A dictionary representing the flattened model.
        """
        obj = {**model.model_dump(), **model.query_index.model_dump()}

        slip_keys = ["query_index"]
        row: ROW = {}
        for key, value in obj.items():
            if key in slip_keys:
                continue
            if isinstance(value, list):
                sub_object = self._flatten_list_field(key, value)
                if sub_object:
                    row.update(sub_object)
            elif isinstance(value, (str, int, float, bool)) or value is None:
                row[key] = value
            else:
                row[key] = str(value)
        if self.fields:
            row = {k: v for k, v in row.items() if k in self.fields}
        return row

    def _flatten_list_field(self, key: str, val: List[Any]) -> ROW:
        """
        Helper method to flatten a list field into a list of strings.

        Args:
            val (List[Any]): The list to flatten.

        Returns:
            List[str]: A list of strings representing the flattened values.
        """
        if not val:
            return {key: []}
        # get the first non-None value
        sample = next((v for v in val if v is not None), None)
        if not sample:
            return {key: []}
        if isinstance(sample, (str, int, float, bool, type(None))):
            return {key: [str(v) for v in val]}
        elif isinstance(sample, dict):
            # TODO: make more schema-driven
            if key.split("_")[-1] in ("terms", "closure", "rollup", "genes"):
                return {
                    f"{key}_id": [v.get('id', str(v)) for v in val],
                    f"{key}_label": [v.get('label', v.get('id', str(v))) for v in val],
                }
            else:
                return {key: [str(v) for v in val]}
        else:
            logger.warning(f"Unexpected type in list: {type(sample)}")
            return {key: [str(v) for v in val]}

