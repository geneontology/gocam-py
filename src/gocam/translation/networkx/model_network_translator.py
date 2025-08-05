from dataclasses import dataclass
from typing import Iterable, Dict, Set, Optional, Any, Union, List
import json

import networkx as nx

from gocam.datamodel import Model, Activity, CausalAssociation
from gocam.translation.networkx.graph_translator import GraphTranslator


@dataclass
class ModelNetworkTranslator(GraphTranslator):
    
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
        # Create query index for efficient lookups
        self.indexer.create_query_index(model)
        
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
        if activity.molecular_function and activity.molecular_function.term:
            attrs[f'{prefix}_molecular_function'] = activity.molecular_function.term
            if activity.molecular_function.evidence:
                # Extract references, evidence codes, and contributors
                references = [e.reference for e in activity.molecular_function.evidence if e.reference]
                evidence_codes = [e.term for e in activity.molecular_function.evidence if e.term]
                contributors = []
                for e in activity.molecular_function.evidence:
                    if e.provenances:
                        for prov in e.provenances:
                            if prov.contributor:
                                contributors.extend(prov.contributor)
                
                if references:
                    attrs[f'{prefix}_molecular_function_has_reference'] = references
                if evidence_codes:
                    attrs[f'{prefix}_molecular_function_assessed_by'] = evidence_codes
                if contributors:
                    attrs[f'{prefix}_molecular_function_contributors'] = contributors
        
        # Add biological process with evidence
        if activity.part_of and activity.part_of.term:
            attrs[f'{prefix}_biological_process'] = activity.part_of.term
            if activity.part_of.evidence:
                references = [e.reference for e in activity.part_of.evidence if e.reference]
                evidence_codes = [e.term for e in activity.part_of.evidence if e.term]
                contributors = []
                for e in activity.part_of.evidence:
                    if e.provenances:
                        for prov in e.provenances:
                            if prov.contributor:
                                contributors.extend(prov.contributor)
                
                if references:
                    attrs[f'{prefix}_biological_process_has_reference'] = references
                if evidence_codes:
                    attrs[f'{prefix}_biological_process_assessed_by'] = evidence_codes
                if contributors:
                    attrs[f'{prefix}_biological_process_contributors'] = contributors
        
        # Add cellular component with evidence
        if activity.occurs_in and activity.occurs_in.term:
            attrs[f'{prefix}_occurs_in'] = activity.occurs_in.term
            if activity.occurs_in.evidence:
                references = [e.reference for e in activity.occurs_in.evidence if e.reference]
                evidence_codes = [e.term for e in activity.occurs_in.evidence if e.term]
                contributors = []
                for e in activity.occurs_in.evidence:
                    if e.provenances:
                        for prov in e.provenances:
                            if prov.contributor:
                                contributors.extend(prov.contributor)
                
                if references:
                    attrs[f'{prefix}_occurs_in_has_reference'] = references
                if evidence_codes:
                    attrs[f'{prefix}_occurs_in_assessed_by'] = evidence_codes
                if contributors:
                    attrs[f'{prefix}_occurs_in_contributors'] = contributors
        
        # Add gene product (enabled_by) with evidence
        if activity.enabled_by and activity.enabled_by.term:
            attrs[f'{prefix}_product'] = activity.enabled_by.term
            if activity.enabled_by.evidence:
                references = [e.reference for e in activity.enabled_by.evidence if e.reference]
                evidence_codes = [e.term for e in activity.enabled_by.evidence if e.term]
                contributors = []
                for e in activity.enabled_by.evidence:
                    if e.provenances:
                        for prov in e.provenances:
                            if prov.contributor:
                                contributors.extend(prov.contributor)
                
                if references:
                    attrs[f'{prefix}_product_has_reference'] = references
                if evidence_codes:
                    attrs[f'{prefix}_product_assessed_by'] = evidence_codes
                if contributors:
                    attrs[f'{prefix}_product_contributors'] = contributors
    
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
