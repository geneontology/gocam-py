from dataclasses import dataclass, field

from gocam.indexing.Indexer import Indexer


@dataclass
class GraphTranslator:
    indexer: Indexer = field(default_factory=lambda: Indexer())