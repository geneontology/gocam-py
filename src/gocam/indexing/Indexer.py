import logging
from dataclasses import dataclass, field
from typing import Optional, List, Collection, Tuple, Any, Set

from oaklib import get_adapter
from oaklib.datamodels.vocabulary import PART_OF, IS_A
from oaklib.interfaces import OboGraphInterface

from gocam.datamodel import Model, Activity, TermAssociation, PublicationObject, TermObject
from gocam.datamodel import QueryIndex
import networkx as nx

logger = logging.getLogger(__name__)

@dataclass
class Indexer:
    """
    Indexes GO-CAM models for querying and analysis.
    
    This class provides methods to:
    1. Index a GO-CAM model by computing statistics and closures
    2. Convert a model to a directed graph
    3. Get term closures for ontology terms
    """
    _go_adapter: OboGraphInterface = None
    subsets: List[str] = field(default_factory=lambda: ["goslim_generic"])

    def go_adapter(self):
        """
        Get or initialize the GO ontology adapter.
        
        Returns:
            An OboGraphInterface implementation for GO
        """
        if not self._go_adapter:
            self._go_adapter = get_adapter("sqlite:obo:go")
        return self._go_adapter
        
    def _get_closures(self, terms: Collection[str]) -> Tuple[List[TermObject], List[TermObject]]:
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
            
        go = self.go_adapter()
        ancs = go.ancestors(list(terms), predicates=[IS_A, PART_OF])
        objs = [
            TermObject(
                id=t,
                label=go.label(t),
            ) for t in terms
        ]
        closure = [
            TermObject(
                id=t,
                label=go.label(t),
            ) for t in ancs if not t.startswith("BFO:")
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


    def create_query_index(self, model: Model, qi: Optional[QueryIndex] = None) -> QueryIndex:
        """
        Create a QueryIndex for the given model.

        :param model:
        :param qi: existing QueryIndex to populate, or None to create a new one
        :return:
        """
        if qi is None:
            qi = QueryIndex()
        go = self.go_adapter()
        model.query_index = qi
        qi = model.query_index
        qi.number_of_activities = len(model.activities or [])
        all_causal_associations = []
        all_refs = set()
        all_mfs = set()
        all_enabled_bys = set()
        all_parts_ofs = set()
        all_occurs_ins = set()
        all_has_inputs = set()
        all_annoton_terms = []

        def term_association_references(term_association: Optional[TermAssociation]) -> Set[str]:
            """
            Extract references from a TermAssociation.

            Returns:
                Set[str]: A set of reference strings extracted from the evidence of the TermAssociation.
            """
            refs = set()
            if term_association and term_association.evidence:
                refs = {e.reference for e in term_association.evidence if e.reference}
            return refs

        for activity in model.activities or []:
            annoton_term_id_parts = []
            if activity.causal_associations:
                all_causal_associations.extend(activity.causal_associations)

            if activity.enabled_by:
                all_refs.update(term_association_references(activity.enabled_by))
                all_enabled_bys.add(activity.enabled_by.term)
                annoton_term_id_parts.append(activity.enabled_by.term)
            
            if activity.molecular_function:
                all_refs.update(term_association_references(activity.molecular_function))
                all_mfs.add(activity.molecular_function.term)
                annoton_term_id_parts.append(activity.molecular_function.term)

            if activity.part_of:
                # Handle both single and list cases
                if isinstance(activity.part_of, list):
                    for ta in activity.part_of:
                        all_refs.update(term_association_references(ta))
                        all_parts_ofs.add(ta.term)
                        annoton_term_id_parts.append(ta.term)
                else:
                    all_refs.update(term_association_references(activity.part_of))
                    all_parts_ofs.add(activity.part_of.term)
                    annoton_term_id_parts.append(activity.part_of.term)
                        
            if activity.occurs_in:
                # Handle both single and list cases
                if isinstance(activity.occurs_in, list):
                    for ta in activity.occurs_in:
                        all_refs.update(term_association_references(ta))
                        all_occurs_ins.add(ta.term)
                        annoton_term_id_parts.append(ta.term)
                else:
                    all_refs.update(term_association_references(activity.occurs_in))
                    all_occurs_ins.add(activity.occurs_in.term)
                    annoton_term_id_parts.append(activity.occurs_in.term)
                
            if activity.has_input:
                for ta in activity.has_input:
                    all_refs.update(term_association_references(ta))
                    all_has_inputs.add(ta.term)
                    annoton_term_id_parts.append(ta.term)

            if activity.enabled_by:
                def _label(x):
                    lbl = go.label(x)
                    if lbl:
                        return lbl
                    else:
                        return x
                annoton_term_label_parts = [_label(x) for x in annoton_term_id_parts]
                annoton_term = TermObject(
                    id="-".join(annoton_term_id_parts),
                    label="; ".join(annoton_term_label_parts),
                )
                all_annoton_terms.append(annoton_term)

        qi.number_of_enabled_by_terms = len(all_enabled_bys)
        qi.number_of_causal_associations = len(all_causal_associations)
        all_refs = list(set(all_refs))
        qi.flattened_references = [
            PublicationObject(id=ref) for ref in all_refs
        ]
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


        subset_terms = set()
        for s in self.subsets:
            subset_terms.update(set(go.subset_members(s)))

        def rollup(terms: Collection[TermObject]) -> List[TermObject]:
            """
            Filter terms to only those in the specified subsets.
            """
            # Use list comprehension instead of set comprehension to avoid hashability issues
            filtered_terms = [t for t in terms if t.id in subset_terms]
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
        qi.model_activity_enabled_by_terms = eb_direct
        qi.model_activity_enabled_by_closure = eb_closure

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





