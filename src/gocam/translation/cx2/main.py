import logging
import re
from typing import Dict, List, Optional, Union

from ndex2.cx2 import CX2Network

from gocam.datamodel import (
    EnabledByProteinComplexAssociation,
    Model,
    MoleculeAssociation,
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


def _remove_species_code_suffix(label: str) -> str:
    for code in SPECIES_CODES:
        label = label.removesuffix(code).strip()
    return label


# Regex from
# https://github.com/ndexbio/ndex-enrichment-rest/wiki/Enrichment-network-structure#via-node-attributes-preferred-method
IQUERY_GENE_SYMBOL_PATTERN = re.compile("(^[A-Z][A-Z0-9-]*$)|(^C[0-9]+orf[0-9]+$)")


def model_to_cx2(gocam: Model) -> list:

    # Internal state
    input_output_nodes: Dict[str, int] = {}
    activity_nodes: Dict[str, int] = {}

    # Internal helper functions that access internal state
    def _get_object_label(object_id: str) -> str:
        object = next((obj for obj in gocam.objects if obj.id == object_id), None)
        return _remove_species_code_suffix(object.label) if object is not None else ""

    def _add_input_output_nodes(
        associations: Optional[Union[MoleculeAssociation, List[MoleculeAssociation]]],
        edge_attributes: dict,
    ) -> None:
        if associations is None:
            return
        if not isinstance(associations, list):
            associations = [associations]
        for association in associations:
            if association.term not in input_output_nodes:
                node_attributes = {
                    "name": _get_object_label(association.term),
                    "represents": association.term,
                    "type": NodeType.MOLECULE.value,
                }

                if association.provenances:
                    node_attributes["provenance"] = [
                        p.contributor for p in association.provenances
                    ]

                input_output_nodes[association.term] = cx2_network.add_node(
                    attributes=node_attributes
                )

            cx2_network.add_edge(
                source=activity_nodes[activity.id],
                target=input_output_nodes[association.term],
                attributes=edge_attributes,
            )

    # Create the CX2 network and set network-level attributes
    cx2_network = CX2Network()
    cx2_network.set_network_attributes(
        {
            "name": gocam.title if gocam.title is not None else gocam.id,
            "represents": gocam.id,
        }
    )

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
            node_type == NodeType.GENE
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
                if IQUERY_GENE_SYMBOL_PATTERN.match(member_name) is None:
                    logger.warning(
                        f"Name for complex member does not match expected pattern: {member_name}"
                    )
                node_attributes["member"].append(member_name)

        if activity.molecular_function:
            node_attributes["molecular_function_id"] = activity.molecular_function.term
            node_attributes["molecular_function_label"] = _get_object_label(
                activity.molecular_function.term
            )

        if activity.occurs_in:
            node_attributes["occurs_in_id"] = activity.occurs_in.term
            node_attributes["occurs_in_label"] = _get_object_label(
                activity.occurs_in.term
            )

        if activity.part_of:
            node_attributes["part_of_id"] = activity.part_of.term
            node_attributes["part_of_label"] = _get_object_label(activity.part_of.term)

        if activity.provenances:
            node_attributes["provenance"] = [
                p.contributor for p in activity.provenances
            ]

        activity_nodes[activity.id] = cx2_network.add_node(attributes=node_attributes)

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
        for association in activity.causal_associations:
            if association.downstream_activity in activity_nodes:
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
                    edge_attributes["evidence"] = [e.term for e in association.evidence]

                if association.provenances:
                    edge_attributes["provenance"] = [
                        p.contributor for p in association.provenances
                    ]

                cx2_network.add_edge(
                    source=activity_nodes[activity.id],
                    target=activity_nodes[association.downstream_activity],
                    attributes=edge_attributes,
                )

    # Set visual properties for the network
    cx2_network.set_visual_properties(VISUAL_PROPERTIES)
    cx2_network.set_opaque_aspect("visualEditorProperties", [VISUAL_EDITOR_PROPERTIES])

    return cx2_network.to_cx2()
