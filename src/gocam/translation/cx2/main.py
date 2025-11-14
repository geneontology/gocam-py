import json
import logging
import re
from functools import cache
from typing import Dict, List, Optional, Union

import networkx as nx
import prefixmaps
from ndex2.cx2 import CX2Network, CX2NetworkXFactory

from gocam.datamodel import (
    EnabledByProteinComplexAssociation,
    EvidenceItem,
    Model,
    MoleculeAssociation,
    TermAssociation,
)
from gocam.translation.cx2.style import (
    RELATIONS,
    VISUAL_EDITOR_PROPERTIES,
    VISUAL_PROPERTIES,
    NodeType,
)
from gocam.utils import remove_species_code_suffix

logger = logging.getLogger(__name__)


# This image gets referenced in the network description, as recommended by NDEx. The process of
# generating this graphic is not fully automated, but it is described here:
# https://github.com/pkalita-lbl/ndex-gocam-legend
LEGEND_GRAPHIC_SRC = "https://geneontology.org/assets/ndex-gocam-legend-v2.png"


@cache
def _get_context():
    return prefixmaps.load_context("go")


def _format_link(url: str, label: str) -> str:
    return f'<a href="{url}" target="_blank" rel="noopener noreferrer">{label}</a>'


# Regex from
# https://github.com/ndexbio/ndex-enrichment-rest/wiki/Enrichment-network-structure#via-node-attributes-preferred-method
IQUERY_GENE_SYMBOL_PATTERN = re.compile("(^[A-Z][A-Z0-9-]*$)|(^C[0-9]+orf[0-9]+$)")


def model_to_cx2(
    gocam: Model, *, validate_iquery_gene_symbol_pattern=True, apply_dot_layout=False
) -> list:
    # Internal state
    input_output_nodes: Dict[str, int] = {}
    activity_nodes_by_activity_id: Dict[str, int] = {}
    activity_nodes_by_enabled_by_id: Dict[str, int] = {}

    go_context = _get_context()
    go_converter = go_context.as_converter()

    # Pre-build object lookup dictionary for better performance
    object_labels = {}
    if gocam.objects:
        for obj in gocam.objects:
            if obj.label:
                object_labels[obj.id] = remove_species_code_suffix(obj.label)
            else:
                object_labels[obj.id] = obj.id

    # Internal helper functions that access internal state
    @cache
    def _get_object_label(object_id: str) -> str:
        return object_labels.get(object_id, object_id)

    @cache
    def _format_curie_link(curie: str) -> str:
        try:
            url = go_converter.expand(curie)
            return _format_link(url, curie)
        except ValueError:
            return curie

    def _format_evidence_list(evidence_list: List[EvidenceItem]) -> str:
        """Format a list of evidence items as an HTML unordered list."""
        if evidence_list is None:
            return ""
        evidence_list_items = []
        for e in evidence_list:
            evidence_item = ""
            if e.reference:
                evidence_item += _format_curie_link(e.reference)
            if e.term:
                term_label = _get_object_label(e.term)
                if evidence_item:
                    evidence_item += f" ({term_label})"
                else:
                    evidence_item += term_label
            if e.with_objects:
                with_objects = ", ".join(_format_curie_link(o) for o in e.with_objects)
                if evidence_item:
                    evidence_item += " "
                evidence_item += f"with/from {with_objects}"
            evidence_list_items.append(f"<li>{evidence_item}</li>")
        return f'<ul style="padding-inline-start: 1rem">{"".join(evidence_list_items)}</ul>'

    def _format_term_association(term_association: TermAssociation) -> str:
        """Format a term association as an HTML link to the term with evidence list."""
        term_id = term_association.term
        term_label = _get_object_label(term_id)
        term_url = go_converter.expand(term_id)
        term_link = _format_link(term_url, f"{term_label} [{term_id}]")
        formatted = f"{term_link}"

        if term_association.evidence:
            evidence_list = _format_evidence_list(term_association.evidence)
            formatted += f"""
<br>
<div style="font-size: smaller; display: block; margin-inline-start: 1rem">
  Evidence:
  {evidence_list}
</div>
"""
        return formatted

    def _add_input_output_nodes(
        associations: Optional[Union[MoleculeAssociation, List[MoleculeAssociation]]],
        edge_attributes: dict,
    ) -> None:
        if associations is None:
            return
        if not isinstance(associations, list):
            associations = [associations]
        for association in associations:
            # Filter proteins at CX2 level (per issue #65)
            # Skip if the term is an INFORMATION_BIOMACROMOLECULE
            # We check if it's already in activity_nodes_by_enabled_by_id as a simple
            # proxy for identifying proteins/gene products
            if (
                association.term in activity_nodes_by_enabled_by_id
                and "has input" in edge_attributes["name"]
            ):
                continue

            if association.term in activity_nodes_by_enabled_by_id:
                target = activity_nodes_by_enabled_by_id[association.term]
            elif association.term in input_output_nodes:
                target = input_output_nodes[association.term]
            else:
                node_attributes = {
                    "name": _get_object_label(association.term),
                    "represents": association.term,
                    "type": NodeType.MOLECULE.value,
                }

                target = cx2_network.add_node(attributes=node_attributes)
                input_output_nodes[association.term] = target

            edge_attributes["Evidence"] = _format_evidence_list(association.evidence)

            cx2_network.add_edge(
                source=activity_nodes_by_activity_id[activity.id],
                target=target,
                attributes=edge_attributes,
            )

    # Create the CX2 network and set network-level attributes
    cx2_network = CX2Network()
    cx2_network.set_network_attributes(
        {
            "@context": json.dumps(go_context.as_dict()),
            "name": gocam.title if gocam.title is not None else gocam.id,
            "prov:wasDerivedFrom": go_converter.expand(gocam.id),
            "description": f'<p><img src="{LEGEND_GRAPHIC_SRC}" style="width: 100%;"/></p>',
        }
    )
    # This gets added separately so we can declare the datatype
    cx2_network.add_network_attribute("labels", [gocam.id], "list_of_string")

    # Add nodes for activities, labeled by the activity's enabled_by object
    for activity in gocam.activities:
        if activity.enabled_by is None:
            continue

        if isinstance(activity.enabled_by, EnabledByProteinComplexAssociation):
            node_type = NodeType.COMPLEX
        else:
            node_type = NodeType.GENE

        node_name = _get_object_label(activity.enabled_by.term)
        if (
            validate_iquery_gene_symbol_pattern
            and node_type == NodeType.GENE
            and IQUERY_GENE_SYMBOL_PATTERN.match(node_name) is None
        ):
            logger.debug(
                f"Name for gene node does not match expected pattern: {node_name}"
            )

        node_attributes = {
            "name": node_name,
            "represents": activity.enabled_by.term,
            "type": node_type.value,
        }

        if node_type == NodeType.COMPLEX and activity.enabled_by.members:
            node_attributes["member"] = []
            for member in activity.enabled_by.members:
                member_name = _get_object_label(member)
                if (
                    validate_iquery_gene_symbol_pattern
                    and IQUERY_GENE_SYMBOL_PATTERN.match(member_name) is None
                ):
                    logger.warning(
                        f"Name for complex member does not match expected pattern: {member_name}"
                    )
                node_attributes["member"].append(member_name)

        node_attributes["Evidence"] = _format_evidence_list(
            activity.enabled_by.evidence if activity.enabled_by else None
        )

        if activity.molecular_function:
            node_attributes["Molecular Function"] = _format_term_association(
                activity.molecular_function
            )

        if activity.occurs_in:
            node_attributes["Occurs In"] = _format_term_association(activity.occurs_in)

        if activity.part_of:
            node_attributes["Part Of"] = _format_term_association(activity.part_of)

        node = cx2_network.add_node(attributes=node_attributes)
        activity_nodes_by_activity_id[activity.id] = node
        activity_nodes_by_enabled_by_id[activity.enabled_by.term] = node

    # Add nodes for input/output molecules and create edges to activity nodes
    for activity in gocam.activities:
        _add_input_output_nodes(
            activity.has_input, {"name": "has input", "represents": "RO:0002233"}
        )
        _add_input_output_nodes(
            activity.has_output, {"name": "has output", "represents": "RO:0002234"}
        )
        _add_input_output_nodes(
            activity.has_primary_input,
            {"name": "has primary input", "represents": "RO:0004009"},
        )
        _add_input_output_nodes(
            activity.has_primary_output,
            {"name": "has primary output", "represents": "RO:0004008"},
        )

    # Add edges for causal associations between activity nodes
    for activity in gocam.activities:
        if activity.causal_associations is None:
            continue

        for association in activity.causal_associations:
            if association.downstream_activity in activity_nodes_by_activity_id:
                relation_style = RELATIONS.get(association.predicate, None)
                if relation_style is None:
                    logger.debug(f"Unknown relation style for {association.predicate}")
                name = (
                    relation_style.label
                    if relation_style is not None
                    else association.predicate
                )
                edge_attributes = {
                    "name": name,
                    "represents": association.predicate,
                }

                if association.evidence:
                    edge_attributes["Evidence"] = _format_evidence_list(
                        association.evidence
                    )

                cx2_network.add_edge(
                    source=activity_nodes_by_activity_id[activity.id],
                    target=activity_nodes_by_activity_id[
                        association.downstream_activity
                    ],
                    attributes=edge_attributes,
                )

    # Set visual properties for the network
    cx2_network.set_visual_properties(VISUAL_PROPERTIES)
    cx2_network.set_opaque_aspect("visualEditorProperties", [VISUAL_EDITOR_PROPERTIES])

    if apply_dot_layout:
        # Convert the CX2 network to a networkx graph
        networkx_graph = CX2NetworkXFactory().get_graph(cx2_network)

        # Our graph, node, and edge attributes confuse the pydot conversion, but we don't need them
        # just for layout
        networkx_graph.graph.clear()
        for node_id in networkx_graph.nodes:
            networkx_graph.nodes[node_id].clear()
        for edge_id in networkx_graph.edges:
            networkx_graph.edges[edge_id].clear()

        # Run graphviz layout on the networkx graph
        layout = nx.nx_pydot.pydot_layout(networkx_graph, prog="dot")

        # These scaling factors are totally heuristic
        x_scale = 2.0
        y_scale = 1.5

        # Get the max x and y values so that we can flip the layout around so its goes top-to-bottom
        # Flipping it left-to-right makes it more similar to the cytoscape-dagre layout used by the
        # pathway widget
        max_x = max([v[0] for v in layout.values()])
        max_y = max([v[1] for v in layout.values()])

        # Stick the computed layout positions back into the CX2 network
        for node_id, position in layout.items():
            cx2_network.get_node(node_id)["x"] = (max_x - position[0]) * x_scale
            cx2_network.get_node(node_id)["y"] = (max_y - position[1]) * y_scale

    return cx2_network.to_cx2()
