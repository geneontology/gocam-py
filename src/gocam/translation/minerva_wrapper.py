import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import DefaultDict, Dict, Iterator, List, Optional, Set, Tuple

import requests
import yaml

from gocam.datamodel import (
    Activity,
    BiologicalProcessAssociation,
    CausalAssociation,
    CellularAnatomicalEntityAssociation,
    EnabledByAssociation,
    EnabledByGeneProductAssociation,
    EnabledByProteinComplexAssociation,
    EvidenceItem,
    Model,
    MolecularFunctionAssociation,
    MoleculeAssociation,
    Object,
    ProvenanceInfo,
)

ENABLED_BY = "RO:0002333"
PART_OF = "BFO:0000050"
HAS_PART = "BFO:0000051"
OCCURS_IN = "BFO:0000066"
HAS_INPUT = "RO:0002233"
HAS_OUTPUT = "RO:0002234"
HAS_PRIMARY_INPUT = "RO:0004009"
HAS_PRIMARY_OUTPUT = "RO:0004008"

logger = logging.getLogger(__name__)


def _normalize_property(prop: str) -> str:
    """
    Normalize a property.

    Sometimes the JSON will use full URIs, sometimes just the local part
    """
    if "/" in prop:
        return prop.split("/")[-1]
    return prop


def _annotations(obj: Dict) -> Dict[str, str]:
    """
    Extract annotations from an object (assumes single-valued).

    Annotations are lists of objects with keys "key" and "value".
    """
    return {_normalize_property(a["key"]): a["value"] for a in obj["annotations"]}


def _annotations_multivalued(obj: Dict) -> Dict[str, List[str]]:
    """
    Extract annotations from an object (assumes multi-valued).

    Annotations are lists of objects with keys "key" and "value".
    """
    anns = defaultdict(list)
    for a in obj["annotations"]:
        anns[a["key"]].append(a["value"])
    return anns


def _setattr_with_warning(obj, attr, value):
    if getattr(obj, attr, None) is not None:
        logger.warning(
            f"Overwriting {attr} for {obj.id if hasattr(obj, 'id') else obj}"
        )
    setattr(obj, attr, value)


MAIN_TYPES = [
    "molecular_function",
    "biological_process",
    "cellular_component",
    "information biomacromolecule",
    "evidence",
    "chemical entity",
    "anatomical entity",
]

COMPLEX_TYPES = [
    "protein-containing complex",
]


@dataclass
class MinervaWrapper:
    """
    An Wrapper over the current GO API which returns "Minerva" JSON objects.

    TODO: Implement a fact counter to ensure all facts are encountered for and nothing dropped on floor
    """

    session: requests.Session = field(default_factory=lambda: requests.Session())
    gocam_index_url: str = "https://go-public.s3.amazonaws.com/files/gocam-models.json"
    gocam_endpoint_base: str = "https://api.geneontology.org/api/go-cam/"

    def models(self) -> Iterator[Model]:
        """Iterator over all GO-CAM models from the index.

        This method fetches the list of all GO-CAM models from the index URL. For each model, the
        Minerva JSON object is fetched and converted to a Model object.

        :return: Iterator over GO-CAM models
        :rtype: Iterator[Model]
        """

        for gocam_id in self.models_ids():
            yield self.fetch_model(gocam_id)

    def models_ids(self) -> Iterator[str]:
        """Iterator over all GO-CAM IDs from the index.

        This method fetches the list of all GO-CAM models from the index URL and returns an
        iterator over the IDs of each model.

        :return: Iterator over GO-CAM IDs
        :rtype: Iterator[str]
        """

        response = self.session.get(self.gocam_index_url)
        response.raise_for_status()
        for model in response.json():
            gocam = model.get("gocam")
            if gocam is None:
                raise ValueError(f"Missing gocam in {model}")
            yield gocam.replace("http://model.geneontology.org/", "")

    def fetch_minerva_object(self, gocam_id: str) -> Dict:
        """Fetch a Minerva JSON object for a given GO-CAM ID.

        :param gocam_id: GO-CAM ID
        :type gocam_id: str
        :return: Minerva JSON object
        :rtype: Dict
        """
        if not gocam_id:
            raise ValueError(f"Missing GO-CAM ID: {gocam_id}")
        local_id = gocam_id.replace("gocam:", "")
        url = f"{self.gocam_endpoint_base}{local_id}"
        logger.info(f"Fetch Minerva JSON from: {url}")
        response = self.session.get(url)
        response.raise_for_status()
        return response.json()

    def fetch_model(self, gocam_id: str) -> Model:
        """Fetch a GO-CAM Model for a given GO-CAM ID.

        :param gocam_id: GO-CAM ID
        :type gocam_id: str
        :return: GO-CAM Model
        :rtype: Model
        """
        minerva_object = self.fetch_minerva_object(gocam_id)
        return self.minerva_object_to_model(minerva_object)

    @staticmethod
    def minerva_object_to_model(obj: Dict) -> Model:
        """Convert a Minerva JSON object to a GO-CAM Model.

        :param obj: Minerva JSON object
        :type obj: Dict
        :return: GO-CAM Model
        :rtype: Model
        """
        id = obj["id"]
        logger.info(f"Processing model id: {id}")
        logger.debug(f"Model: {yaml.dump(obj)}")

        # Bookkeeping variables

        # individual ID to "root" type / category, e.g Evidence, BP
        individual_to_type: Dict[str, Optional[str]] = {}
        individual_to_term: Dict[str, str] = {}
        individual_to_annotations: Dict[str, Dict] = {}
        complex_individuals: Set[str] = set()
        id2obj: Dict[str, Dict] = {}
        activities: List[Activity] = []
        activities_by_mf_id: DefaultDict[str, List[Activity]] = defaultdict(list)
        facts_by_property: DefaultDict[str, List[Dict]] = defaultdict(list)

        def _cls(obj: Dict) -> Optional[str]:
            if obj.get("type", None) == "complement":
                logger.warning(f"Ignoring Complement: {obj}")
                # class expression representing NOT
                return None
            if "id" not in obj:
                raise ValueError(f"No ID for {obj}")
            id = obj["id"]
            id2obj[id] = obj
            return id

        def _evidence_from_fact(fact: Dict) -> List[EvidenceItem]:
            anns_mv = _annotations_multivalued(fact)
            evidence_inst_ids = anns_mv.get("evidence", [])
            evs: List[EvidenceItem] = []
            for evidence_inst_id in evidence_inst_ids:
                evidence_inst_annotations = individual_to_annotations.get(
                    evidence_inst_id, {}
                )
                with_obj: Optional[str] = evidence_inst_annotations.get("with", None)
                if with_obj:
                    with_objs = with_obj.split(" | ")
                else:
                    with_objs = None
                prov = ProvenanceInfo(
                    contributor=evidence_inst_annotations.get("contributor", None),
                    date=evidence_inst_annotations.get("date", None),
                )
                ev = EvidenceItem(
                    term=individual_to_term.get(evidence_inst_id, None),
                    reference=evidence_inst_annotations.get("source", None),
                    with_objects=with_objs,
                    provenances=[prov],
                )
                evs.append(ev)
            return evs

        def _iter_activities_by_fact_subject(
            *,
            fact_property: str,
        ) -> Iterator[Tuple[Activity, str, List[EvidenceItem]]]:
            for fact in facts_by_property.get(fact_property, []):
                s, o = fact["subject"], fact["object"]
                if o not in individual_to_term:
                    logger.warning(f"Missing {o} in {individual_to_term}")
                    continue
                for activity in activities_by_mf_id.get(s, []):
                    evs = _evidence_from_fact(fact)
                    yield activity, individual_to_term[o], evs

        for individual in obj["individuals"]:
            typs = [x["label"] for x in individual.get("root-type", []) if x]
            typ: Optional[str] = None
            for t in typs:
                if t in MAIN_TYPES:
                    typ = t
                    break
            if not typ:
                logger.warning(f"Could not find type for {individual}")
                continue
            individual_to_type[individual["id"]] = typ

            # Check to see if one of the types is a complex type
            for t in typs:
                if t in COMPLEX_TYPES:
                    complex_individuals.add(individual["id"])
                    break

            terms = list(filter(None, (_cls(x) for x in individual.get("type", []))))
            if len(terms) > 1:
                logger.warning(f"Multiple terms for {individual}: {terms}")
            if not terms:
                logger.warning(f"No terms for {individual}")
                continue
            individual_to_term[individual["id"]] = terms[0]
            anns = _annotations(individual)
            individual_to_annotations[individual["id"]] = anns

        for fact in obj["facts"]:
            facts_by_property[fact["property"]].append(fact)

        enabled_by_facts = facts_by_property.get(ENABLED_BY, [])
        if not enabled_by_facts:
            raise ValueError(f"Missing {ENABLED_BY} in {facts_by_property}")
        for fact in enabled_by_facts:
            s, o = fact["subject"], fact["object"]
            if s not in individual_to_term:
                logger.warning(f"Missing {s} in {individual_to_term}")
                continue
            if o not in individual_to_term:
                logger.warning(f"Missing {o} in {individual_to_term}")
                continue
            gene_id = individual_to_term[o]

            evs = _evidence_from_fact(fact)
            enabled_by_association: EnabledByAssociation
            if o in complex_individuals:
                has_part_facts = [
                    fact
                    for fact in facts_by_property.get(HAS_PART, [])
                    if fact["subject"] == o
                ]
                members = [
                    individual_to_term[fact["object"]]
                    for fact in has_part_facts
                    if fact["object"] in individual_to_term
                ]
                enabled_by_association = EnabledByProteinComplexAssociation(
                    term=gene_id, members=members
                )
            else:
                enabled_by_association = EnabledByGeneProductAssociation(term=gene_id)
            activity = Activity(
                id=s,
                enabled_by=enabled_by_association,
                molecular_function=MolecularFunctionAssociation(
                    term=individual_to_term[s], evidence=evs
                ),
            )
            activities.append(activity)
            activities_by_mf_id[s].append(activity)

        for activity, term, evs in _iter_activities_by_fact_subject(
            fact_property=PART_OF
        ):
            association = BiologicalProcessAssociation(term=term, evidence=evs)
            _setattr_with_warning(activity, "part_of", association)

        for activity, term, evs in _iter_activities_by_fact_subject(
            fact_property=OCCURS_IN
        ):
            association = CellularAnatomicalEntityAssociation(term=term, evidence=evs)
            _setattr_with_warning(activity, "occurs_in", association)

        for activity, term, evs in _iter_activities_by_fact_subject(
            fact_property=HAS_INPUT
        ):
            activity.has_input.append(MoleculeAssociation(term=term, evidence=evs))

        for activity, term, evs in _iter_activities_by_fact_subject(
            fact_property=HAS_PRIMARY_INPUT
        ):
            association = MoleculeAssociation(term=term, evidence=evs)
            _setattr_with_warning(activity, "has_primary_input", association)

        for activity, term, evs in _iter_activities_by_fact_subject(
            fact_property=HAS_OUTPUT
        ):
            activity.has_output.append(MoleculeAssociation(term=term, evidence=evs))

        for activity, term, evs in _iter_activities_by_fact_subject(
            fact_property=HAS_PRIMARY_OUTPUT
        ):
            association = MoleculeAssociation(term=term, evidence=evs)
            _setattr_with_warning(activity, "has_primary_output", association)

        for fact_property, facts in facts_by_property.items():
            for fact in facts:
                s, o = fact["subject"], fact["object"]
                subject_activities = activities_by_mf_id.get(s, [])
                object_activities = activities_by_mf_id.get(o, [])

                if not subject_activities or not object_activities:
                    continue
                if individual_to_type.get(s, None) != "molecular_function":
                    continue
                if individual_to_type.get(o, None) != "molecular_function":
                    continue
                if len(subject_activities) > 1:
                    logger.warning(f"Multiple activities for subject: {s}")
                if len(object_activities) > 1:
                    logger.warning(f"Multiple activities for object: {o}")

                subject_activity = subject_activities[0]
                object_activity = object_activities[0]
                evs = _evidence_from_fact(fact)
                rel = CausalAssociation(
                    predicate=fact_property,
                    downstream_activity=object_activity.id,
                    evidence=evs,
                )
                subject_activity.causal_associations.append(rel)

        annotations = _annotations(obj)
        annotations_mv = _annotations_multivalued(obj)
        objs = [Object(id=obj["id"], label=obj["label"]) for obj in id2obj.values()]
        cam = Model(
            id=id,
            title=annotations["title"],
            status=annotations.get("state", None),
            comments=annotations_mv.get("comment", None),
            taxon=annotations.get("in_taxon", None),
            activities=activities,
            objects=objs,
        )
        return cam
