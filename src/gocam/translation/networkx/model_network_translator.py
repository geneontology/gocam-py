from dataclasses import dataclass
from typing import Iterable

import networkx as nx

from gocam.datamodel import Model
from gocam.translation.networkx.graph_translator import GraphTranslator


@dataclass
class ModelNetworkTranslator(GraphTranslator):

    def translate_models(self, models: Iterable[Model]) -> nx.DiGraph:
        for m in models:
            qi = self.indexer.create_query_index(m)
