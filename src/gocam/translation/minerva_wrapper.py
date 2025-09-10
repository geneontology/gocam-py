import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import DefaultDict, Dict, Iterator, List, Optional, Tuple

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
from gocam.vocabulary.taxa import TaxonVocabulary

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
    return {
        _normalize_property(a["key"]): a["value"] for a in obj.get("annotations", [])
    }


def _annotations_multivalued(obj: Dict) -> Dict[str, List[str]]:
    """
    Extract annotations from an object (assumes multi-valued).

    Annotations are lists of objects with keys "key" and "value".
    """
    anns = defaultdict(list)
    for a in obj.get("annotations", []):
        key = _normalize_property(a["key"])
        value = a["value"]
        anns[key].append(value)
    return anns


def _provenance_from_fact(fact: Dict) -> ProvenanceInfo:
    """Produce a ProvenanceInfo object from a fact object."""
    annotations = _annotations(fact)
    annotations_mv = _annotations_multivalued(fact)
    return ProvenanceInfo(
        contributor=annotations_mv.get("contributor"),
        date=annotations.get("date", None),
        provided_by=annotations_mv.get("providedBy"),
    )


def _setattr_with_warning(obj, attr, value):
    if getattr(obj, attr, None) is not None:
        logger.debug(
            f"Overwriting {attr} for {obj.id if hasattr(obj, 'id') else obj}"
        )
    setattr(obj, attr, value)


MOLECULAR_FUNCTION = "GO:0003674"
BIOLOGICAL_PROCESS = "GO:0008150"
CELLULAR_COMPONENT = "GO:0005575"
INFORMATION_BIOMACROMOLECULE = "CHEBI:33695"
PROTEIN_CONTAINING_COMPLEX = "GO:0032991"
EVIDENCE = "ECO:0000000"
CHEMICAL_ENTITY = "CHEBI:24431"
ANATOMICAL_ENTITY = "UBERON:0001062"


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

        # Bookkeeping variables

        # individual ID to "root" type / category, e.g Evidence, BP
        individual_to_root_types: Dict[str, List[str]] = {}
        individual_to_term: Dict[str, str] = {}
        individual_to_annotations: Dict[str, Dict[str, str]] = {}
        individual_to_annotations_multivalued: Dict[str, Dict[str, List[str]]] = {}
        objects_by_id: Dict[str, Dict] = {}
        activities: List[Activity] = []
        activities_by_mf_id: DefaultDict[str, List[Activity]] = defaultdict(list)
        facts_by_property: DefaultDict[str, List[Dict]] = defaultdict(list)

        def _evidence_from_fact(fact: Dict) -> List[EvidenceItem]:
            anns_mv = _annotations_multivalued(fact)
            evidence_inst_ids = anns_mv.get("evidence", [])
            evs: List[EvidenceItem] = []
            for evidence_inst_id in evidence_inst_ids:
                evidence_inst_annotations = individual_to_annotations.get(
                    evidence_inst_id, {}
                )
                evidence_inst_annotations_multivalued = (
                    individual_to_annotations_multivalued.get(evidence_inst_id, {})
                )
                with_obj: Optional[str] = evidence_inst_annotations.get("with", None)
                if with_obj:
                    with_objs = [s.strip() for s in with_obj.split("|")]
                else:
                    with_objs = None
                prov = ProvenanceInfo(
                    contributor=evidence_inst_annotations_multivalued.get(
                        "contributor"
                    ),
                    date=evidence_inst_annotations.get("date", None),
                    provided_by=evidence_inst_annotations_multivalued.get("providedBy"),
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
        ) -> Iterator[Tuple[Activity, str, List[EvidenceItem], ProvenanceInfo]]:
            for fact in facts_by_property.get(fact_property, []):
                subject, object_ = fact["subject"], fact["object"]
                if object_ not in individual_to_term:
                    logger.debug(f"Missing {object_} in {individual_to_term}")
                    continue
                for activity in activities_by_mf_id.get(subject, []):
                    evs = _evidence_from_fact(fact)
                    provenance = _provenance_from_fact(fact)
                    yield activity, object_, evs, provenance

        for individual in obj["individuals"]:
            individual_id = individual["id"]
            root_types = [x["id"] for x in individual.get("root-type", []) if x]
            individual_to_root_types[individual_id] = root_types

            term_id: Optional[str] = None
            for type_ in individual.get("type", []):
                if type_.get("type") == "complement":
                    # class expression representing NOT
                    continue
                type_id = type_.get("id")
                if type_id is None:
                    continue
                objects_by_id[type_id] = type_
                term_id = type_id

            individual_to_term[individual_id] = term_id
            individual_to_annotations[individual_id] = _annotations(individual)
            individual_to_annotations_multivalued[individual_id] = (
                _annotations_multivalued(individual)
            )

        for fact in obj["facts"]:
            facts_by_property[fact["property"]].append(fact)

        enabled_by_facts = facts_by_property.get(ENABLED_BY, [])
        if not enabled_by_facts:
            logger.debug(f"Missing {ENABLED_BY} facts in {facts_by_property}")
        for fact in enabled_by_facts:
            subject, object_ = fact["subject"], fact["object"]
            if subject not in individual_to_term:
                logger.debug(f"Missing {subject} in {individual_to_term}")
                continue
            if object_ not in individual_to_term:
                logger.debug(f"Missing {object_} in {individual_to_term}")
                continue
            gene_id = individual_to_term[object_]
            root_types = individual_to_root_types.get(object_, [])

            evs = _evidence_from_fact(fact)
            prov = _provenance_from_fact(fact)
            enabled_by_association: EnabledByAssociation
            if PROTEIN_CONTAINING_COMPLEX in root_types:
                has_part_facts = [
                    fact
                    for fact in facts_by_property.get(HAS_PART, [])
                    if fact["subject"] == object_
                ]
                members = [
                    individual_to_term[fact["object"]]
                    for fact in has_part_facts
                    if fact["object"] in individual_to_term
                ]
                enabled_by_association = EnabledByProteinComplexAssociation(
                    term=gene_id, members=members, evidence=evs, provenances=[prov]
                )
            elif INFORMATION_BIOMACROMOLECULE in root_types:
                enabled_by_association = EnabledByGeneProductAssociation(
                    term=gene_id, evidence=evs, provenances=[prov]
                )
            else:
                logger.debug(f"Unknown enabled_by type for {object_}; assuming gene product")
                enabled_by_association = EnabledByGeneProductAssociation(
                    term=gene_id, evidence=evs, provenances=[prov]
                )

            activity = Activity(
                id=subject,
                enabled_by=enabled_by_association,
                molecular_function=MolecularFunctionAssociation(
                    term=individual_to_term[subject]
                ),
            )
            activities.append(activity)
            activities_by_mf_id[subject].append(activity)

        for activity, individual, evs, prov in _iter_activities_by_fact_subject(
            fact_property=PART_OF
        ):
            association = BiologicalProcessAssociation(
                term=individual_to_term[individual], evidence=evs, provenances=[prov]
            )
            _setattr_with_warning(activity, "part_of", association)

        for activity, individual, evs, prov in _iter_activities_by_fact_subject(
            fact_property=OCCURS_IN
        ):
            association = CellularAnatomicalEntityAssociation(
                term=individual_to_term[individual], evidence=evs, provenances=[prov]
            )
            _setattr_with_warning(activity, "occurs_in", association)

        for activity, individual, evs, prov in _iter_activities_by_fact_subject(
            fact_property=HAS_INPUT
        ):
            if activity.has_input is None:
                activity.has_input = []
            activity.has_input.append(
                MoleculeAssociation(
                    term=individual_to_term[individual],
                    evidence=evs,
                    provenances=[prov],
                )
            )

        for activity, individual, evs, prov in _iter_activities_by_fact_subject(
            fact_property=HAS_PRIMARY_INPUT
        ):
            association = MoleculeAssociation(
                term=individual_to_term[individual], evidence=evs, provenances=[prov]
            )
            _setattr_with_warning(activity, "has_primary_input", association)

        for activity, individual, evs, prov in _iter_activities_by_fact_subject(
            fact_property=HAS_OUTPUT
        ):
            if activity.has_output is None:
                activity.has_output = []
            activity.has_output.append(
                MoleculeAssociation(
                    term=individual_to_term[individual],
                    evidence=evs,
                    provenances=[prov],
                )
            )

        for activity, individual, evs, prov in _iter_activities_by_fact_subject(
            fact_property=HAS_PRIMARY_OUTPUT
        ):
            association = MoleculeAssociation(
                term=individual_to_term[individual], evidence=evs, provenances=[prov]
            )
            _setattr_with_warning(activity, "has_primary_output", association)

        for fact_property, facts in facts_by_property.items():
            for fact in facts:
                subject, object_ = fact["subject"], fact["object"]
                subject_activities = activities_by_mf_id.get(subject, [])
                object_activities = activities_by_mf_id.get(object_, [])

                if not subject_activities or not object_activities:
                    continue
                if MOLECULAR_FUNCTION not in individual_to_root_types.get(subject, []):
                    continue
                if MOLECULAR_FUNCTION not in individual_to_root_types.get(object_, []):
                    continue
                if len(subject_activities) > 1:
                    logger.debug(f"Multiple activities for subject: {subject}")
                if len(object_activities) > 1:
                    logger.debug(f"Multiple activities for object: {object_}")

                subject_activity = subject_activities[0]
                object_activity = object_activities[0]
                evs = _evidence_from_fact(fact)
                provenance = _provenance_from_fact(fact)
                rel = CausalAssociation(
                    predicate=fact_property,
                    downstream_activity=object_activity.id,
                    evidence=evs,
                    provenances=[provenance],
                )
                if subject_activity.causal_associations is None:
                    subject_activity.causal_associations = []
                subject_activity.causal_associations.append(rel)

        annotations = _annotations(obj)
        annotations_mv = _annotations_multivalued(obj)

        objects: List[Object] = []
        for obj in objects_by_id.values():
            object_ = Object(id=obj["id"])
            if "label" in obj:
                object_.label = obj["label"]
            objects.append(object_)

        provenance = ProvenanceInfo(
            contributor=annotations_mv.get("contributor"),
            date=annotations.get("date", None),
            provided_by=annotations_mv.get("providedBy"),
        )

        # Get all taxa from the annotations
        all_taxa = annotations_mv.get(TaxonVocabulary.TAXON_ANNOTATION_KEY, [])

        # Add legacy taxon key if it exists (backward compatibility)
        legacy_taxon = annotations.get(TaxonVocabulary.LEGACY_TAXON_KEY)
        if legacy_taxon and legacy_taxon not in all_taxa:
            all_taxa.append(legacy_taxon)

        # If no taxa, nothing to do
        if not all_taxa:
            taxon = None
            additional_taxa = []
        # If only one taxon, it's the primary
        elif len(all_taxa) == 1:
            taxon = all_taxa[0]
            additional_taxa = []
        # Multiple taxa - prioritize host taxa as primary
        else:
            # Find host taxa in the list
            host_matches = [t for t in all_taxa if TaxonVocabulary.is_host_taxon(t)]
            if host_matches:
                # Use the first host taxon as primary
                taxon = host_matches[0]
                # All others are additional taxa
                additional_taxa = [t for t in all_taxa if t != taxon]
            else:
                # No host matches, just use the first as primary
                taxon = all_taxa[0]
                additional_taxa = all_taxa[1:]

        # Build model parameters
        model_args = {
            "id": id,
            "title": annotations["title"],
            "status": annotations.get("state", None),
            "comments": annotations_mv.get("comment", None),
            "date_modified": annotations.get("date", None),
            "taxon": taxon,
            "activities": activities,
            "objects": objects,
            "provenances": [provenance],
        }

        # Only add additional_taxa if it has values
        if additional_taxa and len(additional_taxa) > 0:
            model_args["additional_taxa"] = additional_taxa

        cam = Model(**model_args)
        return cam
