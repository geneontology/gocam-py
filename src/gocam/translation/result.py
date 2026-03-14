"""
This module defines classes and enumerations for handling translation results and warnings.

Classes:
    - WarningType: Enumeration of possible translation warning types.
    - TranslationWarning: Represents a warning generated during the translation process.
    - TranslationResult: Holds the result of a translation along with any warnings encountered.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Generic, TypeVar

T = TypeVar("T")


class WarningType(str, Enum):
    """Enumeration of possible translation warning types."""

    ATTRIBUTE_OVERWRITE = "ATTRIBUTE_OVERWRITE"
    INVALID_CIRCULAR_RELATIONSHIP = "INVALID_CIRCULAR_RELATIONSHIP"
    INVALID_MODEL_STATE = "INVALID_MODEL_STATE"
    MISSING_TERM = "MISSING_TERM"
    MULTIPLE_ACTIVITIES = "MULTIPLE_ACTIVITIES"
    NO_ACTIVITIES = "NO_ACTIVITIES"
    NO_ENABLED_BY_FACTS = "NO_ENABLED_BY_FACTS"
    UNHANDLED_FACT = "UNHANDLED_FACT"
    UNKNOWN_ENABLED_BY_TYPE = "UNKNOWN_ENABLED_BY_TYPE"
    UNKNOWN_PROPERTY_INVERSE = "UNKNOWN_PROPERTY_INVERSE"


@dataclass(frozen=True)
class TranslationWarning:
    """Represents a warning generated during the translation process.

    Attributes:
        type: The type or category of the warning.
        message: A descriptive message detailing the warning.
        entity_id: An optional identifier for the entity related to the warning.
    """

    type: WarningType
    message: str
    entity_id: str | None = None


@dataclass(frozen=True)
class TranslationResult(Generic[T]):
    """Holds the result of a translation along with any warnings encountered.

    Attributes:
        result: The translated result object.
        warnings: A list of warnings generated during translation, if any.
    """

    result: T
    warnings: list[TranslationWarning] = field(default_factory=list)
