from dataclasses import dataclass
from typing import Iterable, Dict, Set, Optional, Any, Union, List
import json

import networkx as nx

from gocam.datamodel import Model, Activity, CausalAssociation
from gocam.translation.networkx.graph_translator import GraphTranslator


@dataclass
class ModelNetworkTranslator(GraphTranslator):
    """
    Translates GO-CAM models into gene-to-gene NetworkX DiGraphs.

    This class inherits from GraphTranslator and provides methods to convert one or more 
    GO-CAM Model objects into a NetworkX directed graph where nodes represent gene products 
    and edges represent causal relationships between them.

    Use this class when you need to generate gene-to-gene graphs from GO-CAM models
    with consistent translation interface.
    """
    
    def translate_models(self, models: Iterable[Model]) -> nx.DiGraph:
        """
        Translate multiple GO-CAM models into a single gene-to-gene NetworkX DiGraph.
        
        In the gene-to-gene format:
        - Nodes represent gene products (from enabled_by associations)
        - Edges represent causal relationships between gene products
        - Edge attributes include GO terms (molecular function, biological process, etc.)
        
        Args:
            models: Iterable of GO-CAM Model objects to translate
            
        Returns:
            NetworkX DiGraph where nodes are gene products and edges have GO term properties
        """
        g2g_graph = nx.DiGraph()
        
        for model in models:
            self._add_model_to_graph(model, g2g_graph)
            
        return g2g_graph
    
    def _add_model_to_graph(self, model: Model, graph: nx.DiGraph) -> None:
        """
        Add a single model to the gene-to-gene graph.
        
        Args:
            model: GO-CAM Model to add
            graph: NetworkX DiGraph to add nodes and edges to
        """
        # Build mapping from activity ID to gene product
        activity_to_gene: Dict[str, str] = {}
        
        for activity in model.activities or []:
            if activity.enabled_by and activity.enabled_by.term:
                gene_product = activity.enabled_by.term
                activity_to_gene[activity.id] = gene_product
                
                # Add gene product as node if not already present
                if not graph.has_node(gene_product):
                    node_attrs = self._get_gene_node_attributes(activity, model)
                    graph.add_node(gene_product, **node_attrs)
        
        # Add edges based on causal associations
        for activity in model.activities or []:
            if not activity.causal_associations:
                continue
                
            source_gene = activity_to_gene.get(activity.id)
            if not source_gene:
                continue
                
            for causal_assoc in activity.causal_associations:
                target_gene = activity_to_gene.get(causal_assoc.downstream_activity)
                if not target_gene:
                    continue
                    
                # Create edge with GO term attributes
                edge_attrs = self._get_edge_attributes(
                    activity, 
                    causal_assoc, 
                    model,
                    source_gene,
                    target_gene
                )
                
                # Add or update edge
                if graph.has_edge(source_gene, target_gene):
                    # Merge attributes if edge already exists
                    existing_attrs = graph[source_gene][target_gene]
                    merged_attrs = self._merge_edge_attributes(existing_attrs, edge_attrs)
                    graph[source_gene][target_gene].update(merged_attrs)
                else:
                    graph.add_edge(source_gene, target_gene, **edge_attrs)
    
    def _get_gene_node_attributes(self, activity: Activity, model: Model) -> Dict[str, str]:
        """
        Get node attributes for a gene product.
        
        Args:
            activity: Activity containing the gene product
            model: The GO-CAM model
            
        Returns:
            Dictionary of node attributes
        """
        attrs = {
            'gene_product': activity.enabled_by.term,
            'model_id': model.id
        }
        
        # Add gene product label if available
        if model.objects:
            for obj in model.objects:
                if obj.id == activity.enabled_by.term and obj.label:
                    attrs['label'] = obj.label
                    break
        
        return attrs
    
    def _get_edge_attributes(
        self, 
        source_activity: Activity, 
        causal_assoc: CausalAssociation,
        model: Model,
        source_gene: str,
        target_gene: str
    ) -> Dict[str, Any]:
        """
        Get edge attributes containing GO terms and relationship information.
        
        Args:
            source_activity: Source activity in the causal relationship
            causal_assoc: The causal association
            model: The GO-CAM model
            source_gene: Source gene product ID
            target_gene: Target gene product ID
            
        Returns:
            Dictionary of edge attributes with GO term information
        """
        attrs = {
            'source_gene': source_gene,
            'target_gene': target_gene,
            'model_id': model.id
        }
        
        # Add causal relationship predicate with evidence
        if causal_assoc.predicate:
            attrs['causal_predicate'] = causal_assoc.predicate
            
            # Add evidence for the causal association
            if causal_assoc.evidence:
                references = [e.reference for e in causal_assoc.evidence if e.reference]
                evidence_codes = [e.term for e in causal_assoc.evidence if e.term]
                contributors = []
                for e in causal_assoc.evidence:
                    if e.provenances:
                        for prov in e.provenances:
                            if prov.contributor:
                                contributors.extend(prov.contributor)
                
                if references:
                    attrs['causal_predicate_has_reference'] = references
                if evidence_codes:
                    attrs['causal_predicate_assessed_by'] = evidence_codes
                if contributors:
                    attrs['causal_predicate_contributors'] = contributors
        
        # Add GO terms from source activity
        self._add_activity_go_terms(source_activity, model, attrs, "source_gene")
        
        # Find and add GO terms from target activity
        target_activity = self._find_activity_by_id(causal_assoc.downstream_activity, model)
        if target_activity:
            self._add_activity_go_terms(target_activity, model, attrs, "target_gene")
        
        return attrs
    
    def _find_activity_by_id(self, activity_id: str, model: Model) -> Optional[Activity]:
        """
        Find an activity by its ID in the model.
        
        Args:
            activity_id: The activity ID to search for
            model: The GO-CAM model
            
        Returns:
            Activity object if found, None otherwise
        """
        for activity in model.activities or []:
            if activity.id == activity_id:
                return activity
        return None
    
    def _extract_evidence_data(self, evidence_list):
        """
        Extract references, evidence codes, and contributors from evidence list.
        
        Args:
            evidence_list: List of evidence items
            
        Returns:
            Tuple of (references, evidence_codes, contributors)
        """
        references = [e.reference for e in evidence_list if e.reference]
        evidence_codes = [e.term for e in evidence_list if e.term]
        contributors = []
        for e in evidence_list:
            if e.provenances:
                for prov in e.provenances:
                    if prov.contributor:
                        contributors.extend(prov.contributor)
        return references, evidence_codes, contributors
    
    def _add_term_with_evidence(self, attrs: Dict[str, Any], term_association, prefix: str, term_type: str) -> None:
        """
        Add a term and its evidence to attributes.
        
        Args:
            attrs: Dictionary to add attributes to
            term_association: The term association with evidence
            prefix: Prefix for attribute names ("source_gene" or "target_gene")
            term_type: Type of term (e.g., "molecular_function", "biological_process")
        """
        if term_association and term_association.term:
            attrs[f'{prefix}_{term_type}'] = term_association.term
            if term_association.evidence:
                references, evidence_codes, contributors = self._extract_evidence_data(term_association.evidence)
                
                if references:
                    attrs[f'{prefix}_{term_type}_has_reference'] = references
                if evidence_codes:
                    attrs[f'{prefix}_{term_type}_assessed_by'] = evidence_codes
                if contributors:
                    attrs[f'{prefix}_{term_type}_contributors'] = contributors

    def _add_activity_go_terms(self, activity: Activity, model: Model, attrs: Dict[str, Any], prefix: str) -> None:
        """
        Add GO terms from an activity to the edge attributes with evidence information.
        
        Args:
            activity: The activity to extract GO terms from
            model: The GO-CAM model
            attrs: Dictionary to add attributes to
            prefix: Prefix for attribute names ("source_gene" or "target_gene")
        """
        # Add molecular function with evidence
        self._add_term_with_evidence(attrs, activity.molecular_function, prefix, "molecular_function")
        
        # Add biological process with evidence
        self._add_term_with_evidence(attrs, activity.part_of, prefix, "biological_process")
        
        # Add cellular component with evidence
        self._add_term_with_evidence(attrs, activity.occurs_in, prefix, "occurs_in")
        
        # Add gene product (enabled_by) with evidence
        self._add_term_with_evidence(attrs, activity.enabled_by, prefix, "product")
    
    def _merge_edge_attributes(
        self,
        existing: Dict[str, Union[str, List[str]]],
        new: Dict[str, Union[str, List[str]]]
    ) -> Dict[str, Union[str, List[str]]]:
        """
        Merge edge attributes when multiple causal relationships exist between same genes.
        
        Args:
            existing: Existing edge attributes
            new: New edge attributes to merge
            
        Returns:
            Merged attributes dictionary
        """
        merged = existing.copy()
        
        # For lists of values, we'll concatenate them
        for key, value in new.items():
            if key in merged:
                # Convert to list if not already
                if not isinstance(merged[key], list):
                    merged[key] = [merged[key]]
                if not isinstance(value, list):
                    value = [value]
                
                # Extend the list with new values, avoiding duplicates
                for v in value:
                    if v not in merged[key]:
                        merged[key].append(v)
            else:
                merged[key] = value
                
        return merged
    
    def translate_models_to_json(self, models: Iterable[Model], include_model_info: bool = True, indent: Optional[int] = None) -> str:
        """
        Translate GO-CAM models to gene-to-gene format and return as JSON string.
        
        Args:
            models: Iterable of GO-CAM Model objects to translate
            include_model_info: Whether to include model metadata in the output
            indent: Number of spaces for JSON indentation (None for compact output)
            
        Returns:
            JSON string representation of the gene-to-gene network
        """
        g2g_graph = self.translate_models(models)
        g2g_dict = self._graph_to_dict(g2g_graph, models, include_model_info)
        return json.dumps(g2g_dict, indent=indent)
    
    def _graph_to_dict(self, g2g_graph: nx.DiGraph, models: Iterable[Model], include_model_info: bool) -> Dict:
        """
        Convert NetworkX graph to JSON-serializable dictionary following NetworkX standards.
        
        Args:
            g2g_graph: The gene-to-gene NetworkX DiGraph
            models: The original GO-CAM models (for metadata)
            include_model_info: Whether to include model metadata
            
        Returns:
            Dictionary representation of the gene-to-gene network in NetworkX format
        """
        # Start with NetworkX standard format
        result = {
            "directed": True,
            "multigraph": False,
            "graph": {},
            "nodes": [
                {"id": node, **attrs}
                for node, attrs in g2g_graph.nodes(data=True)
            ],
            "edges": [
                {"source": source, "target": target, **attrs}
                for source, target, attrs in g2g_graph.edges(data=True)
            ]
        }
        
        # Add model metadata to graph attributes following NetworkX standards
        if include_model_info:
            models_list = list(models)
            if len(models_list) == 1:
                # Single model metadata
                model = models_list[0]
                result["graph"]["model_info"] = {
                    "id": model.id,
                    "title": model.title,
                    "taxon": model.taxon,
                    "status": model.status
                }
            else:
                # Multiple models metadata
                result["graph"]["models_info"] = [
                    {
                        "id": model.id,
                        "title": model.title,
                        "taxon": model.taxon,
                        "status": model.status
                    }
                    for model in models_list
                ]
        
        return result
