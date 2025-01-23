import json
import logging
import re
from functools import cache
from typing import Dict, List, Optional, Union

import prefixmaps
from ndex2.cx2 import CX2Network

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

logger = logging.getLogger(__name__)

# Derived from
# https://github.com/geneontology/wc-gocam-viz/blob/6ef1fcaddfef97ece94d04b7c23ac09c33ace168/src/globals/%40noctua.form/data/taxon-dataset.json
# If maintaining this list becomes onerous, consider splitting the label on a space and taking only
# the first part
SPECIES_CODES = [
    "Atal",
    "Btau",
    "Cele",
    "Cfam",
    "Ddis",
    "Dmel",
    "Drer",
    "Ggal",
    "Hsap",
    "Mmus",
    "Pseudomonas",
    "Rnor",
    "Scer",
    "Sjap",
    "Solanaceae",
    "Spom",
    "Sscr",
    "Xenopus",
]

# This graph was produced by the NDEx team based on the style attributes in our CX2 networks. This
# image gets referenced in the network description. It seems a bit fragile to have this static image
# with no process in place to update it if the style changes. But the NDEx folks were fairly
# insistent that we include a legend graphic in the network description.
LEGEND_GRAPHIC_SRC = "https://home.ndexbio.org/img/go-cam_legend_2024108_v2.png"


def _remove_species_code_suffix(label: str) -> str:
    for code in SPECIES_CODES:
        label = label.removesuffix(code).strip()
    return label


@cache
def _get_context():
    return prefixmaps.load_context("go")


def _format_link(url: str, label: str) -> str:
    return f'<a href="{url}" target="_blank" rel="noopener noreferrer">{label}</a>'


# Regex from
# https://github.com/ndexbio/ndex-enrichment-rest/wiki/Enrichment-network-structure#via-node-attributes-preferred-method
IQUERY_GENE_SYMBOL_PATTERN = re.compile("(^[A-Z][A-Z0-9-]*$)|(^C[0-9]+orf[0-9]+$)")


def model_to_cx2(gocam: Model, *, validate_iquery_gene_symbol_pattern=True) -> list:
    # Internal state
    input_output_nodes: Dict[str, int] = {}
    activity_nodes_by_activity_id: Dict[str, int] = {}
    activity_nodes_by_enabled_by_id: Dict[str, int] = {}

    go_context = _get_context()
    go_converter = go_context.as_converter()

    # Internal helper functions that access internal state
    @cache
    def _get_object_label(object_id: str) -> str:
        object = next((obj for obj in gocam.objects if obj.id == object_id), None)
        return _remove_species_code_suffix(object.label) if object is not None else ""

    def _format_evidence_list(evidence_list: List[EvidenceItem]) -> str:
        """Format a list of evidence items as an HTML unordered list."""
        evidence_list_items = []
        for e in evidence_list:
            reference_link = _format_link(go_converter.expand(e.reference), e.reference)
            evidence_item = f"{reference_link} ({_get_object_label(e.term)})"
            if e.with_objects:
                with_objects = ", ".join(
                    _format_link(go_converter.expand(o), o) for o in e.with_objects
                )
                evidence_item += f" with/from {with_objects}"
            evidence_list_items.append(f"<li>{evidence_item}</li>")
        return f'<ul style="padding-inline-start: 1rem">{"".join(evidence_list_items)}</ul>'

    def _format_term_association(term_association: TermAssociation) -> str:
        """Format a term association as an HTML link to the term with evidence list."""
        term_id = term_association.term
        term_label = _get_object_label(term_id)
        term_url = go_converter.expand(term_id)
        term_link = _format_link(term_url, f"{term_label} [{term_id}]")
        evidence_list = _format_evidence_list(term_association.evidence)

        return f"""
{term_link}<br>
<div style="font-size: smaller; display: block; margin-inline-start: 1rem">
  Evidence:
  {evidence_list}
</div>
        """

    def _add_input_output_nodes(
        associations: Optional[Union[MoleculeAssociation, List[MoleculeAssociation]]],
        edge_attributes: dict,
    ) -> None:
        if associations is None:
            return
        if not isinstance(associations, list):
            associations = [associations]
        for association in associations:
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
            "description": f"<p><img src=\"{LEGEND_GRAPHIC_SRC}\" style=\"width: 100%;\"/></p>"
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
            logger.warning(
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
                    logger.warning(
                        f"Unknown relation style for {association.predicate}"
                    )
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

    return cx2_network.to_cx2()
