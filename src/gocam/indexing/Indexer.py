import logging
from collections.abc import Iterable
from functools import cached_property, lru_cache
from typing import Collection, List, Optional, Tuple

import networkx as nx
import pystow
import yaml
from oaklib import get_adapter
from oaklib.datamodels.vocabulary import IS_A, PART_OF

from gocam.datamodel import (
    Activity,
    Association,
    EnabledByGeneProductAssociation,
    EnabledByProteinComplexAssociation,
    Model,
    Object,
    PublicationObject,
    QueryIndex,
    TermObject,
)

logger = logging.getLogger(__name__)

_PYSTOW_MODULE = pystow.module("gocam")
_CURRENT_GROUPS_YAML_URL = "https://current.geneontology.org/metadata/groups.yaml"


def _iter_association_references(
    association: Association | None,
) -> Iterable[str]:
    """
    Extract references from an Association.

    Returns:
        Iterable[str]: Iterable over reference strings extracted from the evidence of the Association.
    """
    if association is not None and association.evidence is not None:
        for evidence in association.evidence:
            if evidence.reference is not None:
                yield evidence.reference


def _iter_provided_bys(
    item: Association | Model | Activity | None,
) -> Iterable[str]:
    """
    Extract all provided_by strings from an item (Association, Model, or Activity) with provenances.

    Returns:
        Iterable[str]: Iterable over provided_by strings extracted from the item's provenances.
    """
    if item is not None and item.provenances is not None:
        for prov in item.provenances:
            if prov.provided_by is not None:
                yield from prov.provided_by


class Indexer:
    """
    Indexes GO-CAM models for querying and analysis.

    This class provides methods to:
    1. Index a GO-CAM model by computing statistics and closures
    2. Convert a model to a directed graph
    3. Get term closures for ontology terms
    """

    def __init__(
        self,
        *,
        subsets: list[str] | None = None,
        go_adapter_descriptor: str = "sqlite:obo:go",
        goc_groups_yaml_path: Optional[str] = None,
        ncbi_taxon_adapter_descriptor: str = "sqlite:obo:ncbitaxon",
    ):
        """
        Initialize the Indexer.

        Args:
            subsets: List of GO subsets to use for rollup. Defaults to ["goslim_generic"].
            go_adapter_descriptor: OAK adapter descriptor for GO ontology. Defaults to "sqlite:obo:go".
            goc_groups_yaml_path: Optional path to a YAML file containing GOC group metadata.
                If not provided, fetches from https://current.geneontology.org/metadata/groups.yaml.
                The YAML file should be a list of group dictionaries, each with an "id" key and other metadata.
            ncbi_taxon_adapter_descriptor: OAK adapter descriptor for NCBI Taxonomy. Defaults to "sqlite:obo:ncbitaxon".
        """
        self._subsets = subsets if subsets is not None else ["goslim_generic"]
        self._go_adapter_descriptor = go_adapter_descriptor
        self._goc_groups_yaml_path = goc_groups_yaml_path
        self._ncbi_taxon_adapter_descriptor = ncbi_taxon_adapter_descriptor

    @cached_property
    def go_adapter(self):
        """
        Get the GO ontology adapter.

        Returns:
            An OboGraphInterface implementation for GO
        """
        return get_adapter(self._go_adapter_descriptor)

    @cached_property
    def ncbi_taxon_adapter(self):
        """
        Get the NCBI Taxonomy ontology adapter.

        Returns:
            An OboGraphInterface implementation for the NCBI Taxonomy database
        """
        return get_adapter(self._ncbi_taxon_adapter_descriptor)

    @cached_property
    def goc_groups(self) -> dict[str, dict]:
        """
        Get metadata about GOC (Gene Ontology Consortium) groups.

        GOC groups are organizational units within the Gene Ontology Consortium,
        each representing a contributing group or database. This property loads
        metadata about these groups from a YAML file, either from a provided path
        or by downloading from the GOC website and caching locally using pystow.

        The returned dictionary maps group IDs (str) to their metadata (dict),
        where each metadata dictionary contains fields such as 'label' and other
        group attributes.

        If the YAML file is unavailable (e.g., the specified path does not exist),
        a FileNotFoundError is raised. If the YAML file is malformed, a
        yaml.YAMLError is raised.

        Caching behavior:
            - The property is cached per Indexer instance via @cached_property.
            - If downloading, the YAML file is cached locally using pystow.

        Returns:
            dict[str, dict]: Dictionary mapping group IDs to metadata dictionaries.

        Raises:
            FileNotFoundError: If the specified goc_groups_yaml_path doesn't exist.
            yaml.YAMLError: If the YAML file is malformed.
        """
        if self._goc_groups_yaml_path is not None:
            path = self._goc_groups_yaml_path
        else:
            path = _PYSTOW_MODULE.ensure(url=_CURRENT_GROUPS_YAML_URL)
        with open(path) as f:
            groups = yaml.safe_load(f)
        return {group["id"]: group for group in groups}

    @cached_property
    def subset_terms(self) -> set[str]:
        """Get all terms that are members of the specified subsets.

        Returns:
            A set of term IDs that are members of the specified subsets.
        """
        subset_terms = set()
        for s in self._subsets:
            subset_terms.update(self.go_adapter.subset_members(s))
        return subset_terms

    @lru_cache(maxsize=None)
    def _go_label(self, id: str) -> str | None:
        """
        Get the label for a GO term.

        Args:
            id: The GO term ID
        """
        return self.go_adapter.label(id)

    @lru_cache(maxsize=None)
    def _ncbi_taxon_label(self, id: str) -> str | None:
        """
        Get the label for an NCBI Taxonomy term.

        Args:
            id: The NCBI Taxonomy term ID
        """
        return self.ncbi_taxon_adapter.label(id)

    def _get_closures(
        self, terms: Collection[str]
    ) -> Tuple[List[TermObject], List[TermObject]]:
        """
        Get direct terms and their transitive closure.

        Args:
            terms: Collection of term IDs to get closures for

        Returns:
            Tuple containing:
            - List of TermObjects for the direct terms
            - List of TermObjects for all ancestors in the closure
        """
        if not terms:
            return [], []

        ancs = self.go_adapter.ancestors(list(terms), predicates=[IS_A, PART_OF])
        objs = [
            TermObject(
                id=t,
                label=self._go_label(t),
            )
            for t in terms
            if t is not None
        ]
        closure = [
            TermObject(
                id=t,
                label=self._go_label(t),
            )
            for t in ancs
            if t is not None and not t.startswith("BFO:")
        ]
        return objs, closure

    def index_model(self, model: Model, reindex=False) -> None:
        """
        Index a GO-CAM model by computing statistics and term closures.

        This method populates the model's query_index with:
        - Basic statistics (number of activities, causal associations)
        - Graph properties (path lengths, strongly connected components)
        - Term closures for molecular functions, biological processes, etc.

        Args:
            model: The GO-CAM model to index
            reindex: Whether to reindex the model if it already has a query_index
        """
        if model.query_index and not reindex:
            return

        if not model.query_index:
            model.query_index = self.create_query_index(model)
        elif reindex:
            self.create_query_index(model, model.query_index)
        else:
            logger.info("Already have index")

    def create_query_index(
        self, model: Model, qi: Optional[QueryIndex] = None
    ) -> QueryIndex:
        """
        Create a QueryIndex for the given model.

        :param model:
        :param qi: existing QueryIndex to populate, or None to create a new one
        :return:
        """
        if qi is None:
            qi = QueryIndex()
        model.query_index = qi
        qi = model.query_index
        qi.number_of_activities = len(model.activities or [])
        all_causal_associations = []
        all_refs = set()
        all_provided_bys: set[str] = set()
        all_mfs = set()
        all_enabled_bys = set()
        all_enabled_by_genes = set()
        all_parts_ofs = set()
        all_occurs_ins = set()
        all_has_inputs = set()
        all_annoton_terms = []
        model_objects_by_id = {obj.id: obj for obj in model.objects or []}

        def _label(x):
            if x in model_objects_by_id and model_objects_by_id[x].label:
                return model_objects_by_id[x].label

            lbl = self._go_label(x)
            if lbl:
                return lbl

            return x

        all_provided_bys.update(_iter_provided_bys(model))

        for activity in model.activities or []:
            all_provided_bys.update(_iter_provided_bys(activity))

            annoton_term_id_parts = []
            if activity.causal_associations:
                all_causal_associations.extend(activity.causal_associations)
                for causal_association in activity.causal_associations:
                    all_refs.update(_iter_association_references(causal_association))
                    all_provided_bys.update(_iter_provided_bys(causal_association))

            if activity.enabled_by:
                all_refs.update(_iter_association_references(activity.enabled_by))
                all_provided_bys.update(_iter_provided_bys(activity.enabled_by))
                all_enabled_bys.add(activity.enabled_by.term)
                if isinstance(activity.enabled_by, EnabledByGeneProductAssociation):
                    all_enabled_by_genes.add(activity.enabled_by.term)
                elif isinstance(
                    activity.enabled_by, EnabledByProteinComplexAssociation
                ):
                    if activity.enabled_by.members:
                        all_enabled_by_genes.update(activity.enabled_by.members)
                annoton_term_id_parts.append(activity.enabled_by.term)

            if activity.molecular_function:
                all_refs.update(
                    _iter_association_references(activity.molecular_function)
                )
                all_provided_bys.update(_iter_provided_bys(activity.molecular_function))
                all_mfs.add(activity.molecular_function.term)
                annoton_term_id_parts.append(activity.molecular_function.term)

            if activity.part_of:
                all_refs.update(_iter_association_references(activity.part_of))
                all_provided_bys.update(_iter_provided_bys(activity.part_of))
                all_parts_ofs.add(activity.part_of.term)
                annoton_term_id_parts.append(activity.part_of.term)

            if activity.occurs_in:
                all_refs.update(_iter_association_references(activity.occurs_in))
                all_provided_bys.update(_iter_provided_bys(activity.occurs_in))
                all_occurs_ins.add(activity.occurs_in.term)
                annoton_term_id_parts.append(activity.occurs_in.term)

            if activity.has_input:
                for ta in activity.has_input:
                    all_refs.update(_iter_association_references(ta))
                    all_provided_bys.update(_iter_provided_bys(ta))
                    all_has_inputs.add(ta.term)
                    annoton_term_id_parts.append(ta.term)

            if activity.enabled_by:
                annoton_term_id_parts = filter(None, annoton_term_id_parts)
                annoton_term_label_parts = [_label(x) for x in annoton_term_id_parts]
                annoton_term = TermObject(
                    id="-".join(annoton_term_id_parts),
                    label="; ".join(annoton_term_label_parts),
                )
                all_annoton_terms.append(annoton_term)

        qi.number_of_enabled_by_terms = len(all_enabled_bys)
        qi.number_of_causal_associations = len(all_causal_associations)
        qi.flattened_references = [PublicationObject(id=ref) for ref in all_refs]
        qi.flattened_provided_by = sorted(
            [
                Object(
                    id=provided_by,
                    label=self.goc_groups.get(provided_by, {}).get(
                        "label", provided_by
                    ),
                )
                for provided_by in all_provided_bys
            ],
            key=lambda x: x.label,
        )
        graph = self.model_to_digraph(model)
        # use nx to find the longest path and all SCCs
        if graph.number_of_nodes() > 0:
            # Find the longest path length
            longest_path = 0
            for node in graph.nodes():
                for other_node in graph.nodes():
                    if node != other_node and nx.has_path(graph, node, other_node):
                        path_length = len(nx.shortest_path(graph, node, other_node)) - 1
                        longest_path = max(longest_path, path_length)
            qi.length_of_longest_causal_association_path = longest_path

            # Create an undirected graph for finding strongly connected components
            # This is because the natural causal graph is a DAG, and we want to
            # find distinct causal subgraphs (weakly connected components)
            undirected_graph = graph.to_undirected()
            connected_components = list(nx.connected_components(undirected_graph))
            qi.number_of_strongly_connected_components = len(connected_components)

        def rollup(terms: Collection[TermObject]) -> List[TermObject]:
            """
            Filter terms to only those in the specified subsets.
            """
            # Use list comprehension instead of set comprehension to avoid hashability issues
            filtered_terms = [t for t in terms if t.id in self.subset_terms]
            # Remove duplicates by id while preserving order
            seen_ids = set()
            result = []
            for t in filtered_terms:
                if t.id not in seen_ids:
                    seen_ids.add(t.id)
                    result.append(t)
            return result

        mf_direct, mf_closure = self._get_closures(all_mfs)
        qi.model_activity_molecular_function_terms = mf_direct
        qi.model_activity_molecular_function_closure = mf_closure
        qi.model_activity_molecular_function_rollup = rollup(mf_closure)

        eb_direct, eb_closure = self._get_closures(all_enabled_bys)
        for term in eb_direct:
            term.label = _label(term.id)
        for term in eb_closure:
            term.label = _label(term.id)
        qi.model_activity_enabled_by_terms = eb_direct
        qi.model_activity_enabled_by_closure = eb_closure

        qi.model_activity_enabled_by_genes = [
            TermObject(id=gene, label=_label(gene))
            for gene in all_enabled_by_genes
            if gene and gene not in ("CHEBI:36080", "CHEBI:33695")
        ]

        parts_ofs_direct, parts_ofs_closure = self._get_closures(all_parts_ofs)
        qi.model_activity_part_of_terms = parts_ofs_direct
        qi.model_activity_part_of_closure = parts_ofs_closure
        qi.model_activity_part_of_rollup = rollup(parts_ofs_closure)

        occurs_in_direct, occurs_in_closure = self._get_closures(all_occurs_ins)
        qi.model_activity_occurs_in_terms = occurs_in_direct
        qi.model_activity_occurs_in_closure = occurs_in_closure
        qi.model_activity_occurs_in_rollup = rollup(occurs_in_closure)

        has_inputs_direct, has_inputs_closure = self._get_closures(all_has_inputs)
        qi.model_activity_has_input_terms = has_inputs_direct
        qi.model_activity_has_input_closure = has_inputs_closure
        qi.model_activity_has_input_rollup = rollup(has_inputs_closure)

        qi.annoton_terms = all_annoton_terms

        if model.taxon:
            qi.taxon_label = self._ncbi_taxon_label(model.taxon)

        return qi

    def model_to_digraph(self, model: Model) -> nx.DiGraph:
        """
        Convert a model to a directed graph where nodes are activities
        and edges represent causal relationships between activities.

        Args:
            model: The GO-CAM model to convert

        Returns:
            A directed graph (DiGraph) where nodes are activity IDs and edges represent
            causal relationships from source to target activities
        """
        g = nx.DiGraph()
        for a in model.activities or []:
            if a.causal_associations:
                for ca in a.causal_associations:
                    if ca.downstream_activity:
                        g.add_edge(ca.downstream_activity, a.id)
        return g
