from dataclasses import dataclass
from typing import Iterable, Dict, Set, Optional

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
        qi = self.indexer.create_query_index(model)
        
        # Build mapping from activity ID to gene product
        activity_to_gene: Dict[str, str] = {}
        gene_to_activities: Dict[str, Set[str]] = {}
        
        for activity in model.activities or []:
            if activity.enabled_by and activity.enabled_by.term:
                gene_product = activity.enabled_by.term
                activity_to_gene[activity.id] = gene_product
                
                if gene_product not in gene_to_activities:
                    gene_to_activities[gene_product] = set()
                gene_to_activities[gene_product].add(activity.id)
                
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
    ) -> Dict[str, str]:
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
        
        # Add causal relationship predicate
        if causal_assoc.predicate:
            attrs['causal_predicate'] = causal_assoc.predicate
            
            # Add predicate label if available
            if model.objects:
                for obj in model.objects:
                    if obj.id == causal_assoc.predicate and obj.label:
                        attrs['causal_predicate_label'] = obj.label
                        break
        
        # Add molecular function from source activity
        if source_activity.molecular_function and source_activity.molecular_function.term:
            attrs['source_molecular_function'] = source_activity.molecular_function.term
            
            # Add MF label if available
            if model.objects:
                for obj in model.objects:
                    if obj.id == source_activity.molecular_function.term and obj.label:
                        attrs['source_molecular_function_label'] = obj.label
                        break
        
        # Add biological process from source activity
        if source_activity.part_of and source_activity.part_of.term:
            attrs['source_biological_process'] = source_activity.part_of.term
            
            # Add BP label if available  
            if model.objects:
                for obj in model.objects:
                    if obj.id == source_activity.part_of.term and obj.label:
                        attrs['source_biological_process_label'] = obj.label
                        break
        
        # Add cellular component from source activity
        if source_activity.occurs_in and source_activity.occurs_in.term:
            attrs['source_occurs_in'] = source_activity.occurs_in.term
            
            # Add CC label if available
            if model.objects:
                for obj in model.objects:
                    if obj.id == source_activity.occurs_in.term and obj.label:
                        attrs['source_occurs_in_label'] = obj.label
                        break
        
        return attrs
    
    def _merge_edge_attributes(self, existing: Dict, new: Dict) -> Dict:
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
