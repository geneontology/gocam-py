import logging
from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import Any

import requests

from gocam.datamodel import (
    Activity,
    BiologicalProcessAssociation,
    CausalAssociation,
    CellTypeAssociation,
    CellularAnatomicalEntityAssociation,
    EnabledByAssociation,
    EnabledByGeneProductAssociation,
    EnabledByProteinComplexAssociation,
    EvidenceItem,
    Model,
    ModelStateEnum,
    MolecularFunctionAssociation,
    MoleculeAssociation,
    MoleculeNode,
    Object,
    PhaseAssociation,
    ProteinComplexMemberAssociation,
    ProvenanceInfo,
)
from gocam.translation.result import TranslationResult, TranslationWarning, WarningType
from gocam.vocabulary import Relation, TaxonVocabulary

MOLECULAR_FUNCTION = "GO:0003674"
BIOLOGICAL_PROCESS = "GO:0008150"
CELLULAR_COMPONENT = "GO:0005575"
INFORMATION_BIOMACROMOLECULE = "CHEBI:33695"
PROTEIN_CONTAINING_COMPLEX = "GO:0032991"
EVIDENCE = "ECO:0000000"
CHEMICAL_ENTITY = "CHEBI:24431"
ANATOMICAL_ENTITY = "UBERON:0001062"

# These are properties which occur in facts that represent associations between molecular functions
# and molecules, where the subject is the molecular function and the object is the molecule. These
# facts can be directly translated into MoleculeAssociation objects on the relevant Activity.
MOLECULAR_ASSOCIATION_PROPERTIES = (
    Relation.HAS_INPUT,
    Relation.HAS_PRIMARY_INPUT,
    Relation.HAS_OUTPUT,
    Relation.HAS_PRIMARY_OUTPUT,
    Relation.HAS_SMALL_MOLECULE_ACTIVATOR,
    Relation.HAS_SMALL_MOLECULE_INHIBITOR,
)

# These keys of this dict are properties which occur in facts that represent associations between
# molecular functions and molecules, where the subject is the molecule and the object is the
# molecular function. The values of the dict are the corresponding inverse properties. Because the
# data model is oriented around activities (molecular functions), if we encounter a fact with one of
# the key properties, we will use the corresponding inverse property to create a MoleculeAssociation
# on the relevant Activity.
MOLECULAR_ASSOCIATION_INVERSE_PROPERTIES = {
    Relation.INPUT_OF: Relation.HAS_INPUT,
    Relation.OUTPUT_OF: Relation.HAS_OUTPUT,
    Relation.IS_SMALL_MOLECULE_ACTIVATOR_OF: Relation.HAS_SMALL_MOLECULE_ACTIVATOR,
    Relation.IS_SMALL_MOLECULE_INHIBITOR_OF: Relation.HAS_SMALL_MOLECULE_INHIBITOR,
}

logger = logging.getLogger(__name__)


@dataclass
class MinervaTranslator:
    """
    Translates a Minerva JSON object to a GO-CAM Model.
    """

    minerva_obj: dict
    individual_to_root_types: dict[str, list[str]] = field(default_factory=dict)
    individual_to_term: dict[str, str] = field(default_factory=dict)
    individual_to_annotations: dict[str, dict[str, str]] = field(default_factory=dict)
    individual_to_annotations_multivalued: dict[str, dict[str, list[str]]] = field(
        default_factory=dict
    )
    objects_by_id: dict[str, dict] = field(default_factory=dict)
    activities: list[Activity] = field(default_factory=list)
    activities_by_mf_id: defaultdict[str, list[Activity]] = field(
        default_factory=lambda: defaultdict(list)
    )
    molecule_nodes_by_id: dict[str, MoleculeNode] = field(default_factory=dict)
    facts_by_property: defaultdict[str, list[dict]] = field(
        default_factory=lambda: defaultdict(list)
    )
    facts_by_subject_property: defaultdict[tuple[str, str], list[dict]] = field(
        default_factory=lambda: defaultdict(list)
    )
    processed_facts: set[tuple[str, str, str]] = field(default_factory=set)
    translation_warnings: set[TranslationWarning] = field(default_factory=set)

    def translate(self) -> TranslationResult[Model]:
        """
        Execute the translation process.

        Returns:
            Object containing GO-CAM Model and any translation warnings
        """
        self._parse_individuals()
        self._parse_facts()
        self._create_activities()
        self._process_biological_process_associations()
        self._process_occurs_in_associations()
        self._process_happens_during_associations()
        self._process_molecule_associations()
        self._process_causal_associations()

        return self._build_result()

    def _normalize_property(self, prop: str) -> str:
        """Normalize a property URI."""
        if "/" in prop:
            return prop.split("/")[-1]
        return prop

    def _annotations(self, obj: dict) -> dict[str, str]:
        """Extract single-valued annotations."""
        return {
            self._normalize_property(a["key"]): a["value"]
            for a in obj.get("annotations", [])
        }

    def _annotations_multivalued(self, obj: dict) -> dict[str, list[str]]:
        """Extract multi-valued annotations."""
        anns = defaultdict(list)
        for a in obj.get("annotations", []):
            key = self._normalize_property(a["key"])
            value = a["value"]
            anns[key].append(value)
        return anns

    def _setattr_with_warning(self, obj: Activity, attr: str, value: Any):
        """Set an attribute on an Activity, with a warning if it already exists."""
        if getattr(obj, attr, None) is not None:
            self.translation_warnings.add(
                TranslationWarning(
                    type=WarningType.ATTRIBUTE_OVERWRITE,
                    message=f"Overwriting {attr} for {obj.id}",
                    entity_id=obj.id,
                )
            )
        setattr(obj, attr, value)

    def _process_fact(self, fact: dict) -> tuple[list[EvidenceItem], ProvenanceInfo]:
        """Process a fact to extract evidence and provenance information.

        This method also marks the fact as processed by adding it to the processed_facts set, which
        is used to track which facts have been handled during translation and identify any that
        were not.
        """
        self.processed_facts.add((fact["subject"], fact["property"], fact["object"]))

        annotations = self._annotations(fact)
        annotations_mv = self._annotations_multivalued(fact)

        evidence_inst_ids = annotations_mv.get("evidence", [])
        evs: list[EvidenceItem] = []
        for evidence_inst_id in evidence_inst_ids:
            evidence_inst_annotations = self.individual_to_annotations.get(
                evidence_inst_id, {}
            )
            evidence_inst_annotations_multivalued = (
                self.individual_to_annotations_multivalued.get(evidence_inst_id, {})
            )
            with_obj: str | None = evidence_inst_annotations.get("with", None)
            if with_obj:
                with_objs = [s.strip() for s in with_obj.split("|")]
            else:
                with_objs = None
            prov = ProvenanceInfo(
                contributor=evidence_inst_annotations_multivalued.get("contributor"),
                date=evidence_inst_annotations.get("date", None),
                provided_by=evidence_inst_annotations_multivalued.get("providedBy"),
            )
            ev = EvidenceItem(
                term=self.individual_to_term.get(evidence_inst_id, None),
                reference=evidence_inst_annotations.get("source", None),
                with_objects=with_objs,
                provenances=[prov],
            )
            evs.append(ev)

        prov = ProvenanceInfo(
            contributor=annotations_mv.get("contributor"),
            date=annotations.get("date", None),
            provided_by=annotations_mv.get("providedBy"),
        )
        return evs, prov

    def _iter_activities_by_fact_subject(
        self,
        *,
        fact_property: str,
    ) -> Iterator[tuple[Activity, str, list[EvidenceItem], ProvenanceInfo]]:
        """Iterate over activities that are the subject of a fact with the given property."""
        for fact in self.facts_by_property.get(fact_property, []):
            subject, object_ = fact["subject"], fact["object"]
            if object_ not in self.individual_to_term:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.MISSING_TERM,
                        message=f"Missing term for object {object_} in fact",
                        entity_id=object_,
                    )
                )
                continue
            for activity in self.activities_by_mf_id.get(subject, []):
                evs, provenance = self._process_fact(fact)
                yield activity, object_, evs, provenance

    def _add_molecule_node(self, individual_id: str) -> MoleculeNode | None:
        """Add a molecule node for a given individual ID, if it doesn't already exist."""
        if individual_id in self.molecule_nodes_by_id:
            return self.molecule_nodes_by_id[individual_id]

        if individual_id not in self.individual_to_term:
            self.translation_warnings.add(
                TranslationWarning(
                    type=WarningType.MISSING_TERM,
                    message=f"Missing term for individual {individual_id} when creating molecule node",
                    entity_id=individual_id,
                )
            )
            return None

        molecule_node = MoleculeNode(
            id=individual_id,
            term=self.individual_to_term[individual_id],
        )

        for located_in_fact in self.facts_by_subject_property.get(
            (individual_id, Relation.LOCATED_IN), []
        ):
            evs, prov = self._process_fact(located_in_fact)
            if molecule_node.located_in is not None:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.ATTRIBUTE_OVERWRITE,
                        message=f"Overwriting located_in for molecule {individual_id}",
                        entity_id=individual_id,
                    )
                )
            molecule_node.located_in = CellularAnatomicalEntityAssociation(
                term=self.individual_to_term.get(located_in_fact["object"], None),
                evidence=evs,
                provenances=[prov],
            )

        self.molecule_nodes_by_id[individual_id] = molecule_node
        return molecule_node

    def _add_molecule_association(
        self, fact: dict, *, is_inverse: bool = False
    ) -> None:
        """Add a molecule association to the relevant activity based on a fact."""
        if not is_inverse:
            mf_individual_id = fact["subject"]
            molecule_individual_id = fact["object"]
            predicate = fact["property"]
        else:
            mf_individual_id = fact["object"]
            molecule_individual_id = fact["subject"]
            predicate = MOLECULAR_ASSOCIATION_INVERSE_PROPERTIES.get(
                fact["property"], None
            )
            if predicate is None:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.UNKNOWN_PROPERTY_INVERSE,
                        message=f"Fact property '{fact['property']}' does not have a known inverse",
                        entity_id=fact["property"],
                    )
                )
                return

        molecule_node = self._add_molecule_node(molecule_individual_id)
        if molecule_node is None:
            return

        activities = self.activities_by_mf_id.get(mf_individual_id)
        if not activities:
            self.translation_warnings.add(
                TranslationWarning(
                    type=WarningType.NO_ACTIVITIES,
                    message=f"No activities found for molecular function individual {mf_individual_id}",
                    entity_id=mf_individual_id,
                )
            )
            return
        evs, prov = self._process_fact(fact)
        for activity in activities:
            association = MoleculeAssociation(
                molecule=molecule_node.id,
                predicate=predicate,
                evidence=evs,
                provenances=[prov],
            )
            if activity.molecular_associations is None:
                activity.molecular_associations = []
            activity.molecular_associations.append(association)

    def _parse_individuals(self):
        """Parse individuals from the Minerva object."""
        for individual in self.minerva_obj["individuals"]:
            individual_id = individual["id"]
            root_types = [x["id"] for x in individual.get("root-type", []) if x]
            self.individual_to_root_types[individual_id] = root_types

            term_id: str | None = None
            for type_ in individual.get("type", []):
                if type_.get("type") == "complement":
                    # class expression representing NOT
                    continue
                type_id = type_.get("id")
                if type_id is None:
                    continue
                self.objects_by_id[type_id] = type_
                term_id = type_id

            if term_id is None:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.MISSING_TERM,
                        message=f"Missing term for individual {individual_id}",
                        entity_id=individual_id,
                    )
                )
                continue

            self.individual_to_term[individual_id] = term_id
            self.individual_to_annotations[individual_id] = self._annotations(
                individual
            )
            self.individual_to_annotations_multivalued[individual_id] = (
                self._annotations_multivalued(individual)
            )

    def _parse_facts(self):
        """Parse facts from the Minerva object."""
        for fact in self.minerva_obj["facts"]:
            self.facts_by_property[fact["property"]].append(fact)
            self.facts_by_subject_property[(fact["subject"], fact["property"])].append(
                fact
            )

    def _create_activities(self):
        """Create Activity objects from enabled_by facts."""
        enabled_by_facts = self.facts_by_property.get(Relation.ENABLED_BY, [])
        if not enabled_by_facts:
            self.translation_warnings.add(
                TranslationWarning(
                    type=WarningType.NO_ENABLED_BY_FACTS,
                    message=f"No enabled_by ({Relation.ENABLED_BY}) facts found",
                    entity_id=self.minerva_obj["id"],
                )
            )
        for enabled_by_fact in enabled_by_facts:
            subject, object_ = enabled_by_fact["subject"], enabled_by_fact["object"]
            if subject not in self.individual_to_term:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.MISSING_TERM,
                        message=f"Missing term for subject {subject} in fact",
                        entity_id=subject,
                    )
                )
                continue
            if object_ not in self.individual_to_term:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.MISSING_TERM,
                        message=f"Missing term for object {object_} in fact",
                        entity_id=object_,
                    )
                )
                continue
            gene_id = self.individual_to_term[object_]
            root_types = self.individual_to_root_types.get(object_, [])

            evs, prov = self._process_fact(enabled_by_fact)
            enabled_by_association: EnabledByAssociation
            if PROTEIN_CONTAINING_COMPLEX in root_types:
                member_associations: list[ProteinComplexMemberAssociation] = []
                has_part_facts = self.facts_by_subject_property.get(
                    (object_, Relation.HAS_PART), []
                )
                for has_part_fact in has_part_facts:
                    member_object = has_part_fact["object"]
                    member_term = self.individual_to_term.get(member_object)
                    if member_term is None:
                        self.translation_warnings.add(
                            TranslationWarning(
                                type=WarningType.MISSING_TERM,
                                message=f"Missing term for object {member_object} in has_part fact",
                                entity_id=member_object,
                            )
                        )
                        continue
                    member_evs, member_prov = self._process_fact(has_part_fact)
                    member_associations.append(
                        ProteinComplexMemberAssociation(
                            term=member_term,
                            evidence=member_evs,
                            provenances=[member_prov],
                        )
                    )
                enabled_by_association = EnabledByProteinComplexAssociation(
                    term=gene_id,
                    members=member_associations,
                    evidence=evs,
                    provenances=[prov],
                )
            elif INFORMATION_BIOMACROMOLECULE in root_types:
                enabled_by_association = EnabledByGeneProductAssociation(
                    term=gene_id, evidence=evs, provenances=[prov]
                )
            else:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.UNKNOWN_ENABLED_BY_TYPE,
                        message=f"Unknown enabled_by type for {object_}; assuming gene product",
                        entity_id=object_,
                    )
                )
                enabled_by_association = EnabledByGeneProductAssociation(
                    term=gene_id, evidence=evs, provenances=[prov]
                )

            activity = Activity(
                id=subject,
                enabled_by=enabled_by_association,
                molecular_function=MolecularFunctionAssociation(
                    term=self.individual_to_term[subject]
                ),
            )
            self.activities.append(activity)
            self.activities_by_mf_id[subject].append(activity)

    def _process_biological_process_associations(self):
        """Process part_of facts to build biological process associations."""
        for (
            activity,
            object_,
            evs,
            prov,
        ) in self._iter_activities_by_fact_subject(fact_property=Relation.PART_OF):
            association = BiologicalProcessAssociation(
                term=self.individual_to_term[object_], evidence=evs, provenances=[prov]
            )

            happens_during_facts = self.facts_by_subject_property.get(
                (object_, Relation.HAPPENS_DURING), []
            )
            for happens_during_fact in happens_during_facts:
                evs, prov = self._process_fact(happens_during_fact)
                if association.happens_during is not None:
                    self.translation_warnings.add(
                        TranslationWarning(
                            type=WarningType.ATTRIBUTE_OVERWRITE,
                            message=f"Overwriting part_of.happens_during for {activity.id}",
                            entity_id=activity.id,
                        )
                    )
                association.happens_during = PhaseAssociation(
                    term=self.individual_to_term.get(
                        happens_during_fact["object"], None
                    ),
                    evidence=evs,
                    provenances=[prov],
                )

            part_of_facts = self.facts_by_subject_property.get(
                (object_, Relation.PART_OF), []
            )
            for part_of_fact in part_of_facts:
                evs, prov = self._process_fact(part_of_fact)
                if association.part_of is not None:
                    self.translation_warnings.add(
                        TranslationWarning(
                            type=WarningType.ATTRIBUTE_OVERWRITE,
                            message=f"Overwriting part_of.part_of for {activity.id}",
                            entity_id=activity.id,
                        )
                    )
                association.part_of = BiologicalProcessAssociation(
                    term=self.individual_to_term.get(part_of_fact["object"], None),
                    evidence=evs,
                    provenances=[prov],
                )

            self._setattr_with_warning(activity, "part_of", association)

    def _process_occurs_in_associations(self):
        """Process occurs_in facts to build cellular anatomical entity associations."""
        for (
            activity,
            object_,
            evs,
            prov,
        ) in self._iter_activities_by_fact_subject(fact_property=Relation.OCCURS_IN):
            association = CellularAnatomicalEntityAssociation(
                term=self.individual_to_term[object_], evidence=evs, provenances=[prov]
            )

            part_of_facts = self.facts_by_subject_property.get(
                (object_, Relation.PART_OF), []
            )
            for part_of_fact in part_of_facts:
                evs, prov = self._process_fact(part_of_fact)
                if association.part_of is not None:
                    self.translation_warnings.add(
                        TranslationWarning(
                            type=WarningType.ATTRIBUTE_OVERWRITE,
                            message=f"Overwriting occurs_in.part_of for {activity.id}",
                            entity_id=activity.id,
                        )
                    )
                association.part_of = CellTypeAssociation(
                    term=self.individual_to_term.get(part_of_fact["object"], None),
                    evidence=evs,
                    provenances=[prov],
                )

            self._setattr_with_warning(activity, "occurs_in", association)

    def _process_happens_during_associations(self):
        """Process happens_during facts to build phase associations."""
        for (
            activity,
            object_,
            evs,
            prov,
        ) in self._iter_activities_by_fact_subject(
            fact_property=Relation.HAPPENS_DURING
        ):
            association = PhaseAssociation(
                term=self.individual_to_term[object_], evidence=evs, provenances=[prov]
            )
            self._setattr_with_warning(activity, "happens_during", association)

    def _process_molecule_associations(self):
        """Process facts representing associations between molecular functions and molecules."""
        for fact_property, facts in self.facts_by_property.items():
            if fact_property in MOLECULAR_ASSOCIATION_PROPERTIES:
                for fact in facts:
                    self._add_molecule_association(fact)
            elif fact_property in MOLECULAR_ASSOCIATION_INVERSE_PROPERTIES:
                for fact in facts:
                    self._add_molecule_association(fact, is_inverse=True)

    def _process_causal_associations(self):
        """Process facts representing causal associations between activities."""
        for fact_property, facts in self.facts_by_property.items():
            for fact in facts:
                subject, object_ = fact["subject"], fact["object"]
                subject_activities = self.activities_by_mf_id.get(subject, [])
                object_activities = self.activities_by_mf_id.get(object_, [])

                if not subject_activities or not object_activities:
                    continue
                if MOLECULAR_FUNCTION not in self.individual_to_root_types.get(
                    subject, []
                ):
                    continue
                if MOLECULAR_FUNCTION not in self.individual_to_root_types.get(
                    object_, []
                ):
                    continue
                if len(subject_activities) > 1:
                    self.translation_warnings.add(
                        TranslationWarning(
                            type=WarningType.MULTIPLE_ACTIVITIES,
                            message=f"Multiple activities for subject {subject}",
                            entity_id=subject,
                        )
                    )
                if len(object_activities) > 1:
                    self.translation_warnings.add(
                        TranslationWarning(
                            type=WarningType.MULTIPLE_ACTIVITIES,
                            message=f"Multiple activities for object {object_}",
                            entity_id=object_,
                        )
                    )

                subject_activity = subject_activities[0]
                object_activity = object_activities[0]
                evs, provenance = self._process_fact(fact)
                rel = CausalAssociation(
                    predicate=fact_property,
                    downstream_activity=object_activity.id,
                    evidence=evs,
                    provenances=[provenance],
                )
                if subject_activity.causal_associations is None:
                    subject_activity.causal_associations = []
                subject_activity.causal_associations.append(rel)

    def _build_result(self) -> TranslationResult[Model]:
        """Build the final TranslationResult object."""
        annotations = self._annotations(self.minerva_obj)
        annotations_mv = self._annotations_multivalued(self.minerva_obj)

        objects: list[Object] = []
        for obj in self.objects_by_id.values():
            object_ = Object(id=obj["id"])
            if "label" in obj:
                object_.label = obj["label"]
            objects.append(object_)

        provenance = ProvenanceInfo(
            contributor=annotations_mv.get("contributor"),
            date=annotations.get("date", None),
            provided_by=annotations_mv.get("providedBy"),
        )

        all_taxa = annotations_mv.get(TaxonVocabulary.TAXON_ANNOTATION_KEY, [])
        legacy_taxon = annotations.get(TaxonVocabulary.LEGACY_TAXON_KEY)
        if legacy_taxon and legacy_taxon not in all_taxa:
            all_taxa.append(legacy_taxon)

        if not all_taxa:
            taxon = None
            additional_taxa = []
        elif len(all_taxa) == 1:
            taxon = all_taxa[0]
            additional_taxa = []
        else:
            host_matches = [t for t in all_taxa if TaxonVocabulary.is_host_taxon(t)]
            if host_matches:
                taxon = host_matches[0]
                additional_taxa = [t for t in all_taxa if t != taxon]
            else:
                taxon = all_taxa[0]
                additional_taxa = all_taxa[1:]

        num_processed_facts = len(self.processed_facts)
        num_facts = len(self.minerva_obj["facts"])
        if num_processed_facts < num_facts:
            for fact in self.minerva_obj["facts"]:
                fact_tuple = (fact["subject"], fact["property"], fact["object"])
                if fact_tuple not in self.processed_facts:
                    self.translation_warnings.add(
                        TranslationWarning(
                            type=WarningType.UNHANDLED_FACT,
                            message=f"Unhandled fact: {fact_tuple}",
                            entity_id=self.minerva_obj["id"],
                        )
                    )

        cam = Model(
            id=self.minerva_obj["id"],
            title=annotations["title"],
            comments=annotations_mv.get("comment", None),
            date_modified=annotations.get("date", None),
            taxon=taxon,
            activities=self.activities,
            molecules=list(self.molecule_nodes_by_id.values()),
            objects=objects,
            provenances=[provenance],
        )

        state_annotation = annotations.get("state", None)
        if state_annotation:
            try:
                cam.status = ModelStateEnum(state_annotation)
            except ValueError:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.INVALID_MODEL_STATE,
                        message=f"Invalid model state {state_annotation}",
                        entity_id=self.minerva_obj["id"],
                    )
                )

        if additional_taxa and len(additional_taxa) > 0:
            cam.additional_taxa = additional_taxa

        return TranslationResult(result=cam, warnings=list(self.translation_warnings))


@dataclass
class MinervaWrapper:
    """
    A Wrapper over the current GO API which returns "Minerva" JSON objects.
    """

    session: requests.Session = field(default_factory=lambda: requests.Session())
    gocam_index_url: str = "https://go-public.s3.amazonaws.com/files/gocam-models.json"
    gocam_endpoint_base: str = "https://api.geneontology.org/api/go-cam/"

    def models(self) -> Iterator[Model]:
        """Iterator over all GO-CAM models from the index.

        This method fetches the list of all GO-CAM models from the index URL. For each model, the
        Minerva JSON object is fetched and converted to a Model object.

        Yields:
            GO-CAM Model
        """

        for gocam_id in self.models_ids():
            yield self.fetch_model(gocam_id)

    def models_ids(self) -> Iterator[str]:
        """Iterator over all GO-CAM IDs from the index.

        This method fetches the list of all GO-CAM models from the index URL and returns an
        iterator over the IDs of each model.

        Yields:
            GO-CAM ID
        """

        response = self.session.get(self.gocam_index_url)
        response.raise_for_status()
        for model in response.json():
            gocam = model.get("gocam")
            if gocam is None:
                raise ValueError(f"Missing gocam in {model}")
            yield gocam.replace("http://model.geneontology.org/", "")

    def fetch_minerva_object(self, gocam_id: str) -> dict:
        """Fetch a Minerva JSON object for a given GO-CAM ID.

        Args:
            gocam_id: GO-CAM ID

        Returns:
            Minerva JSON object
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

        Args:
            gocam_id: GO-CAM ID

        Returns:
            GO-CAM Model
        """
        minerva_object = self.fetch_minerva_object(gocam_id)
        return self.minerva_object_to_model(minerva_object)

    @staticmethod
    def translate(minerva_obj: dict) -> TranslationResult[Model]:
        """Convert a Minerva JSON object to a GO-CAM Model.

        Args:
            minerva_obj: Minerva JSON object

        Returns:
            Object containing GO-CAM Model and any translation warnings
        """
        translator = MinervaTranslator(minerva_obj)
        return translator.translate()

    @staticmethod
    def minerva_object_to_model(obj: dict) -> Model:
        """Convert a Minerva JSON object to a GO-CAM Model.

        Args:
            obj: Minerva JSON object

        Returns:
            GO-CAM Model
        """
        result = MinervaWrapper.translate(obj)
        return result.result
