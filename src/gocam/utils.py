# Derived from:
# https://github.com/geneontology/web-components/blob/5d87e593121eafe6ac4690fa4591f88aa5a03fd8/packages/web-components/src/globals/%40noctua.form/data/taxon-dataset.json
from collections import defaultdict

import networkx as nx

from gocam.datamodel import Activity, Model, MoleculeAssociation

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


def remove_species_code_suffix(label: str) -> str:
    """Remove known species codes from the end of a label.

    If a label ends with one of the known species codes, remove it and any trailing whitespace.
    Otherwise, return the label unchanged.

    :param label: The label to process.
    :return: The processed label.
    """
    for code in SPECIES_CODES:
        label = label.removesuffix(code).strip()
    return label


def all_activity_inputs(activity: Activity) -> list[MoleculeAssociation]:
    """Get all inputs for an activity, including primary input if present.

    Args:
        activity: The activity to get inputs for.

    Returns:
        List of all molecule associations that are inputs to the activity.
    """
    inputs: list[MoleculeAssociation] = list(activity.has_input or [])
    if activity.has_primary_input:
        inputs.append(activity.has_primary_input)
    return inputs


def all_activity_outputs(activity: Activity) -> list[MoleculeAssociation]:
    """Get all outputs for an activity, including primary output if present.

    Args:
        activity: The activity to get outputs for.

    Returns:
        List of all molecule associations that are outputs of the activity.
    """
    outputs: list[MoleculeAssociation] = list(activity.has_output or [])
    if activity.has_primary_output:
        outputs.append(activity.has_primary_output)
    return outputs


def model_to_digraph(model: Model) -> nx.DiGraph:
    """Convert a Model to a NetworkX directed graph.

    Nodes are added for activities with an 'enabled_by' property. Edges are added based on:
    - Causal associations between activities.
    - Outputs of one activity serving as inputs to another activity.

    Args:
        model: The model to convert.

    Returns:
        A directed graph representation of the model.
    """
    graph = nx.DiGraph()

    activities_by_input: dict[str, list[str]] = defaultdict(list)
    for activity in model.activities or []:
        for input_ in all_activity_inputs(activity):
            if input_.molecule:
                activities_by_input[input_.molecule].append(activity.id)

    for activity in model.activities or []:
        if activity.enabled_by is None:
            continue

        graph.add_node(activity.id)

        for causal_assoc in activity.causal_associations or []:
            downstream_activity_id = causal_assoc.downstream_activity
            if downstream_activity_id:
                graph.add_edge(activity.id, downstream_activity_id)

        for output in all_activity_outputs(activity):
            if output.molecule:
                for downstream_activity_id in activities_by_input.get(
                    output.molecule, []
                ):
                    if downstream_activity_id != activity.id:
                        graph.add_edge(activity.id, downstream_activity_id)

    return graph
