from ndex2.cx2 import CX2Network

from gocam.datamodel import EnabledByProteinComplexAssociation, Model
from gocam.translation.cx2.style import (
    RELATIONS,
    VISUAL_EDITOR_PROPERTIES,
    VISUAL_PROPERTIES,
)


def model_to_cx2(gocam: Model) -> list:

    def _get_object_label(object_id: str) -> str:
        object = next((obj for obj in gocam.objects if obj.id == object_id), None)
        return object.label if object is not None else ""

    cx2_network = CX2Network()
    cx2_network.set_network_attributes(
        {
            "name": gocam.title if gocam.title is not None else gocam.id,
            "represents": gocam.id,
        }
    )

    activity_nodes = {}
    for activity in gocam.activities:
        if activity.enabled_by is None:
            continue

        if isinstance(activity.enabled_by, EnabledByProteinComplexAssociation):
            node_type = "complex"
        else:
            node_type = "gene"

        node_attributes = {
            "name": _get_object_label(activity.enabled_by.term),
            "represents": activity.enabled_by.term,
            "type": node_type,
        }

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

    for activity in gocam.activities:
        for association in activity.causal_associations:
            if association.downstream_activity in activity_nodes:
                relation_style = RELATIONS.get(association.predicate, None)
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

    cx2_network.set_visual_properties(VISUAL_PROPERTIES)
    cx2_network.set_opaque_aspect("visualEditorProperties", [VISUAL_EDITOR_PROPERTIES])

    return cx2_network.to_cx2()
