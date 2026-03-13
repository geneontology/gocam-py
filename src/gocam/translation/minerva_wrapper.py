import logging
from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import Any

import requests

from gocam.datamodel import (
    Activity,
    Association,
    BiologicalProcessAssociation,
    CausalAssociation,
    CellTypeAssociation,
    CellularAnatomicalEntityAssociation,
    EnabledByAssociation,
    EnabledByGeneProductAssociation,
    EnabledByProteinComplexAssociation,
    EvidenceItem,
    GrossAnatomyAssociation,
    Model,
    ModelStateEnum,
    MolecularFunctionAssociation,
    MoleculeAssociation,
    MoleculeNode,
    Object,
    PhaseAssociation,
    ProteinComplexMemberAssociation,
    ProvenanceInfo,
    TermAssociation,
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
class MinervaView:
    """
    Indexed, read-only view of raw Minerva JSON data.
    """

    raw_json: dict
    _facts_by_property: defaultdict[str, list[dict]] = field(
        init=False, default_factory=lambda: defaultdict(list)
    )
    _facts_by_subject_property: defaultdict[tuple[str, str], list[dict]] = field(
        init=False, default_factory=lambda: defaultdict(list)
    )
    _facts_by_object_property: defaultdict[tuple[str, str], list[dict]] = field(
        init=False, default_factory=lambda: defaultdict(list)
    )
    _individual_id_to_term: dict[str, str] = field(init=False, default_factory=dict)
    _individual_id_to_root_types: dict[str, list[str]] = field(
        init=False, default_factory=dict
    )
    _individual_id_to_annotations: dict[str, dict[str, str]] = field(
        init=False, default_factory=dict
    )
    _individual_id_to_annotations_multivalued: dict[str, dict[str, list[str]]] = field(
        init=False, default_factory=dict
    )
    _objects_by_id: dict[str, dict] = field(init=False, default_factory=dict)

    def __post_init__(self):
        self._facts_by_property = defaultdict(list)
        self._facts_by_subject_property = defaultdict(list)
        self._facts_by_object_property = defaultdict(list)

        for individual in self.raw_json.get("individuals", []):
            ind_id = individual["id"]
            self._individual_id_to_root_types[ind_id] = [
                x["id"] for x in individual.get("root-type", []) if x
            ]
            self._individual_id_to_annotations[ind_id] = self._extract_annotations(
                individual
            )
            self._individual_id_to_annotations_multivalued[ind_id] = (
                self._extract_annotations_multivalued(individual)
            )

            for type_ in individual.get("type", []):
                if type_.get("type") == "complement":
                    continue
                type_id = type_.get("id")
                if type_id:
                    self._objects_by_id[type_id] = type_
                    self._individual_id_to_term[ind_id] = type_id

        for fact in self.raw_json.get("facts", []):
            prop = fact["property"]
            subj = fact["subject"]
            obj = fact["object"]
            self._facts_by_property[prop].append(fact)
            self._facts_by_subject_property[(subj, prop)].append(fact)
            self._facts_by_object_property[(obj, prop)].append(fact)

    def _normalize_property(self, prop: str) -> str:
        """Normalize a property URI."""
        if "/" in prop:
            return prop.split("/")[-1]
        return prop

    def _extract_annotations(self, obj: dict) -> dict[str, str]:
        """Extract single-valued annotations."""
        return {
            self._normalize_property(a["key"]): a["value"]
            for a in obj.get("annotations", [])
        }

    def _extract_annotations_multivalued(self, obj: dict) -> dict[str, list[str]]:
        """Extract multi-valued annotations."""
        anns = defaultdict(list)
        for a in obj.get("annotations", []):
            key = self._normalize_property(a["key"])
            value = a["value"]
            anns[key].append(value)
        return anns

    def get_term(self, individual_id: str) -> str | None:
        """Get the term (GO, CHEBI, ECO, etc.) associated with an individual."""
        return self._individual_id_to_term.get(individual_id)

    def get_root_types(self, individual_id: str) -> list[str]:
        """Get the root types associated with an individual."""
        return self._individual_id_to_root_types.get(individual_id, [])

    def is_type(self, individual_id: str, type_uri: str) -> bool:
        """Check if an individual has a specific root type."""
        return type_uri in self.get_root_types(individual_id)

    def get_annotations(self, id_or_dict: str | dict) -> dict[str, str]:
        """Get single-valued annotations for an individual (by ID) or a fact (by dict)."""
        if isinstance(id_or_dict, str):
            return self._individual_id_to_annotations.get(id_or_dict, {})
        return self._extract_annotations(id_or_dict)

    def get_annotations_multivalued(
        self, id_or_dict: str | dict
    ) -> dict[str, list[str]]:
        """Get multivalued annotations for an individual (by ID) or a fact (by dict)."""
        if isinstance(id_or_dict, str):
            return self._individual_id_to_annotations_multivalued.get(id_or_dict, {})
        return self._extract_annotations_multivalued(id_or_dict)

    def get_facts(
        self,
        subject: str | None = None,
        object: str | None = None,
        property: str | None = None,
    ) -> list[dict]:
        """Query facts by subject, object, and/or predicate."""
        if subject and property:
            return self._facts_by_subject_property.get((subject, property), [])
        if object and property:
            return self._facts_by_object_property.get((object, property), [])
        if property:
            return self._facts_by_property.get(property, [])
        # Fallback to iteration for other filters
        return [
            f
            for f in self.raw_json.get("facts", [])
            if (not subject or f["subject"] == subject)
            and (not object or f["object"] == object)
        ]

    def all_objects(self) -> list[dict]:
        """Get all object definitions."""
        return list(self._objects_by_id.values())

    def all_facts(self) -> list[dict]:
        """Get all facts."""
        return self.raw_json.get("facts", [])


@dataclass
class MinervaTranslator:
    """
    Translates a Minerva JSON object to a GO-CAM Model.
    """

    minerva_obj: dict
    view: MinervaView = field(init=False)
    activities: list[Activity] = field(default_factory=list)
    activities_by_mf_id: defaultdict[str, list[Activity]] = field(
        default_factory=lambda: defaultdict(list)
    )
    molecule_nodes_by_id: dict[str, MoleculeNode] = field(default_factory=dict)
    processed_facts: set[tuple[str, str, str]] = field(default_factory=set)
    translation_warnings: set[TranslationWarning] = field(default_factory=set)

    def __post_init__(self):
        self.view = MinervaView(self.minerva_obj)

    def translate(self) -> TranslationResult[Model]:
        """
        Execute the translation process.

        Returns:
            Object containing GO-CAM Model and any translation warnings
        """
        self._create_activities()
        self._process_biological_process_associations()
        self._process_occurs_in_associations()
        self._process_happens_during_associations()
        self._process_molecule_associations()
        self._process_causal_associations()

        return self._build_result()

    def _set_activity_attr(self, activity: Activity, attr: str, value: Any):
        """Set an attribute on an Activity, with a warning if it already exists."""
        if getattr(activity, attr, None) is not None:
            self.translation_warnings.add(
                TranslationWarning(
                    type=WarningType.ATTRIBUTE_OVERWRITE,
                    message=f"Overwriting {attr} for activity {activity.id}",
                    entity_id=activity.id,
                )
            )
        setattr(activity, attr, value)

    def _set_term_association_attr(
        self, association: TermAssociation, attr: str, value: Any
    ):
        """Set an attribute on an TermAssociation, with a warning if it already exists."""
        if getattr(association, attr, None) is not None:
            self.translation_warnings.add(
                TranslationWarning(
                    type=WarningType.ATTRIBUTE_OVERWRITE,
                    message=f"Overwriting {attr} for association with term {association.term}",
                    entity_id=association.term,
                )
            )
        setattr(association, attr, value)

    def _process_fact(self, fact: dict) -> tuple[list[EvidenceItem], ProvenanceInfo]:
        """Process a fact to extract evidence and provenance information.

        This method also marks the fact as processed by adding it to the processed_facts set, which
        is used to track which facts have been handled during translation and identify any that
        were not.
        """
        self.processed_facts.add((fact["subject"], fact["property"], fact["object"]))
        return self._extract_evidence_and_provenance(fact)

    def _extract_evidence_and_provenance(
        self, id_or_dict: str | dict
    ) -> tuple[list[EvidenceItem], ProvenanceInfo]:
        """Extract evidence and provenance for an entity or fact."""
        annotations = self.view.get_annotations(id_or_dict)
        annotations_mv = self.view.get_annotations_multivalued(id_or_dict)

        evidence_individual_ids = annotations_mv.get("evidence", [])
        evs: list[EvidenceItem] = []
        for evidence_individual_id in evidence_individual_ids:
            evidence_annotations = self.view.get_annotations(evidence_individual_id)
            evidence_annotations_mv = self.view.get_annotations_multivalued(
                evidence_individual_id
            )

            with_obj: str | None = evidence_annotations.get("with", None)
            if with_obj:
                with_objs = [s.strip() for s in with_obj.split("|")]
            else:
                with_objs = None

            prov = ProvenanceInfo(
                contributor=evidence_annotations_mv.get("contributor"),
                date=evidence_annotations.get("date"),
                provided_by=evidence_annotations_mv.get("providedBy"),
            )
            ev = EvidenceItem(
                term=self.view.get_term(evidence_individual_id),
                reference=evidence_annotations.get("source"),
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

    def _create_activities(self):
        """Create Activity objects from enabled_by facts."""
        enabled_by_facts = self.view.get_facts(property=Relation.ENABLED_BY)
        if not enabled_by_facts:
            self.translation_warnings.add(
                TranslationWarning(
                    type=WarningType.NO_ENABLED_BY_FACTS,
                    message=f"No enabled_by ({Relation.ENABLED_BY}) facts found",
                    entity_id=self.minerva_obj["id"],
                )
            )
        for enabled_by_fact in enabled_by_facts:
            enabled_by_subject = enabled_by_fact["subject"]
            enabled_by_object = enabled_by_fact["object"]
            subject_term = self.view.get_term(enabled_by_subject)
            object_term = self.view.get_term(enabled_by_object)

            if subject_term is None:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.MISSING_TERM,
                        message=f"Missing term for subject {enabled_by_subject} in fact",
                        entity_id=enabled_by_subject,
                    )
                )
                continue
            if object_term is None:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.MISSING_TERM,
                        message=f"Missing term for object {enabled_by_object} in fact",
                        entity_id=enabled_by_object,
                    )
                )
                continue

            evs, prov = self._process_fact(enabled_by_fact)
            enabled_by_association = self._build_enabled_by_association(
                enabled_by_object, object_term, evs, prov
            )

            activity = Activity(
                id=enabled_by_subject,
                enabled_by=enabled_by_association,
                molecular_function=MolecularFunctionAssociation(term=subject_term),
            )
            self.activities.append(activity)
            self.activities_by_mf_id[enabled_by_subject].append(activity)

    def _build_enabled_by_association(
        self,
        individual_id: str,
        term: str,
        evidence: list[EvidenceItem],
        prov: ProvenanceInfo,
    ) -> EnabledByAssociation:
        """Build an EnabledByAssociation (GeneProduct or ProteinComplex)."""
        if self.view.is_type(individual_id, PROTEIN_CONTAINING_COMPLEX):
            return self._build_protein_complex_association(
                individual_id, term, evidence, prov
            )

        if not self.view.is_type(individual_id, INFORMATION_BIOMACROMOLECULE):
            self.translation_warnings.add(
                TranslationWarning(
                    type=WarningType.UNKNOWN_ENABLED_BY_TYPE,
                    message=f"Unknown enabled_by type for {individual_id}; assuming gene product",
                    entity_id=individual_id,
                )
            )

        return EnabledByGeneProductAssociation(
            term=term, evidence=evidence, provenances=[prov]
        )

    def _build_protein_complex_association(
        self,
        individual_id: str,
        term: str,
        evidence: list[EvidenceItem],
        prov: ProvenanceInfo,
    ) -> EnabledByProteinComplexAssociation:
        """Build an EnabledByProteinComplexAssociation by gathering member facts."""
        member_associations: list[ProteinComplexMemberAssociation] = []
        for has_part_fact in self.view.get_facts(
            subject=individual_id, property=Relation.HAS_PART
        ):
            member_id = has_part_fact["object"]
            member_term = self.view.get_term(member_id)
            if member_term is None:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.MISSING_TERM,
                        message=f"Missing term for object {member_id} in has_part fact",
                        entity_id=member_id,
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
        return EnabledByProteinComplexAssociation(
            term=term,
            members=member_associations,
            evidence=evidence,
            provenances=[prov],
        )

    def _add_molecule_node(self, individual_id: str) -> MoleculeNode | None:
        """Add a molecule node for a given individual ID, if it doesn't already exist."""
        if individual_id in self.molecule_nodes_by_id:
            return self.molecule_nodes_by_id[individual_id]

        term = self.view.get_term(individual_id)
        if term is None:
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
            term=term,
        )

        for located_in_fact in self.view.get_facts(
            subject=individual_id, property=Relation.LOCATED_IN
        ):
            if molecule_node.located_in is not None:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.ATTRIBUTE_OVERWRITE,
                        message=f"Overwriting located_in for molecule {individual_id}",
                        entity_id=individual_id,
                    )
                )
            molecule_node.located_in = (
                self._build_cellular_anatomical_entity_association(located_in_fact)
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

    def _build_biological_process_association(
        self, fact: dict, visited: frozenset[str] | None = None
    ) -> BiologicalProcessAssociation | None:
        """Recursively build a BiologicalProcessAssociation."""
        individual_id = fact["object"]

        if visited is None:
            visited = frozenset()
        if individual_id in visited:
            self.translation_warnings.add(
                TranslationWarning(
                    type=WarningType.INVALID_CIRCULAR_RELATIONSHIP,
                    message=f"Circular relationship detected for {individual_id}",
                    entity_id=individual_id,
                )
            )
            return None
        next_visited = visited.union({individual_id})

        term = self.view.get_term(individual_id)
        evidence, provenance = self._process_fact(fact)
        association = BiologicalProcessAssociation(
            term=term,
            evidence=evidence,
            provenances=[provenance],
        )

        for happens_during_fact in self.view.get_facts(
            subject=individual_id, property=Relation.HAPPENS_DURING
        ):
            phase_association = self._build_phase_association(happens_during_fact)
            self._set_term_association_attr(
                association, "happens_during", phase_association
            )

        for part_of_fact in self.view.get_facts(
            subject=individual_id, property=Relation.PART_OF
        ):
            bp_association = self._build_biological_process_association(
                part_of_fact, visited=next_visited
            )
            if bp_association is not None:
                self._set_term_association_attr(association, "part_of", bp_association)

        return association

    def _build_phase_association(self, fact: dict) -> PhaseAssociation:
        """Build a PhaseAssociation."""
        individual_id = fact["object"]
        evidence, provenance = self._process_fact(fact)
        term = self.view.get_term(individual_id)
        return PhaseAssociation(
            term=term,
            evidence=evidence,
            provenances=[provenance],
        )

    def _build_cellular_anatomical_entity_association(
        self, fact: dict
    ) -> CellularAnatomicalEntityAssociation:
        """Build a CellularAnatomicalEntityAssociation."""
        individual_id = fact["object"]
        evidence, provenance = self._process_fact(fact)
        term = self.view.get_term(individual_id)
        association = CellularAnatomicalEntityAssociation(
            term=term,
            evidence=evidence,
            provenances=[provenance],
        )

        for part_of_fact in self.view.get_facts(
            subject=individual_id, property=Relation.PART_OF
        ):
            cell_type_association = self._build_cell_type_association(part_of_fact)
            self._set_term_association_attr(
                association, "part_of", cell_type_association
            )

        return association

    def _build_cell_type_association(self, fact: dict) -> CellTypeAssociation:
        """Build a CellTypeAssociation."""
        individual_id = fact["object"]
        evidence, provenance = self._process_fact(fact)
        term = self.view.get_term(individual_id)
        association = CellTypeAssociation(
            term=term,
            evidence=evidence,
            provenances=[provenance],
        )

        for part_of_fact in self.view.get_facts(
            subject=individual_id, property=Relation.PART_OF
        ):
            gross_anatomy_association = self._build_gross_anatomy_association(
                part_of_fact
            )
            self._set_term_association_attr(
                association, "part_of", gross_anatomy_association
            )

        return association

    def _build_gross_anatomy_association(
        self, fact: dict, visited: frozenset[str] | None = None
    ) -> GrossAnatomyAssociation | None:
        """Recursively build a GrossAnatomyAssociation."""
        individual_id = fact["object"]

        if visited is None:
            visited = frozenset()
        if individual_id in visited:
            self.translation_warnings.add(
                TranslationWarning(
                    type=WarningType.INVALID_CIRCULAR_RELATIONSHIP,
                    message=f"Circular relationship detected for {individual_id}",
                    entity_id=individual_id,
                )
            )
            return None
        next_visited = visited.union({individual_id})

        evidence, provenance = self._process_fact(fact)
        term = self.view.get_term(individual_id)
        association = GrossAnatomyAssociation(
            term=term,
            evidence=evidence,
            provenances=[provenance],
        )

        for part_of_fact in self.view.get_facts(
            subject=individual_id, property=Relation.PART_OF
        ):
            gross_anatomy_association = self._build_gross_anatomy_association(
                part_of_fact,
                visited=next_visited,
            )
            if gross_anatomy_association is not None:
                self._set_term_association_attr(
                    association, "part_of", gross_anatomy_association
                )

        return association

    def _process_biological_process_associations(self):
        """Process part_of facts to build biological process associations."""
        for fact in self.view.get_facts(property=Relation.PART_OF):
            part_of_subject = fact["subject"]
            part_of_object = fact["object"]

            activities = self.activities_by_mf_id.get(part_of_subject)
            if not activities:
                continue

            if self.view.get_term(part_of_object) is None:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.MISSING_TERM,
                        message=f"Missing term for object {part_of_object} in part_of fact",
                        entity_id=part_of_object,
                    )
                )
                continue

            association = self._build_biological_process_association(fact)

            for activity in activities:
                self._set_activity_attr(activity, "part_of", association)

    def _process_occurs_in_associations(self):
        """Process occurs_in facts to build cellular anatomical entity associations."""
        for fact in self.view.get_facts(property=Relation.OCCURS_IN):
            occurs_in_subject = fact["subject"]
            occurs_in_object = fact["object"]

            activities = self.activities_by_mf_id.get(occurs_in_subject)
            if not activities:
                continue

            if self.view.get_term(occurs_in_object) is None:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.MISSING_TERM,
                        message=f"Missing term for object {occurs_in_object} in occurs_in fact",
                        entity_id=occurs_in_object,
                    )
                )
                continue

            association = self._build_cellular_anatomical_entity_association(fact)

            for activity in activities:
                self._set_activity_attr(activity, "occurs_in", association)

    def _process_happens_during_associations(self):
        """Process happens_during facts to build phase associations."""
        for fact in self.view.get_facts(property=Relation.HAPPENS_DURING):
            happens_during_subject = fact["subject"]
            happens_during_object = fact["object"]

            activities = self.activities_by_mf_id.get(happens_during_subject)
            if not activities:
                continue

            if self.view.get_term(happens_during_object) is None:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.MISSING_TERM,
                        message=f"Missing term for object {happens_during_object} in happens_during fact",
                        entity_id=happens_during_object,
                    )
                )
                continue

            association = self._build_phase_association(fact)

            for activity in activities:
                self._set_activity_attr(activity, "happens_during", association)

    def _process_molecule_associations(self):
        """Process facts representing associations between molecular functions and molecules."""
        for fact in self.view.all_facts():
            fact_property = fact["property"]
            if fact_property in MOLECULAR_ASSOCIATION_PROPERTIES:
                self._add_molecule_association(fact)
            elif fact_property in MOLECULAR_ASSOCIATION_INVERSE_PROPERTIES:
                self._add_molecule_association(fact, is_inverse=True)

    def _process_causal_associations(self):
        """Process facts representing causal associations between activities."""
        for fact in self.view.all_facts():
            fact_property = fact["property"]
            fact_subject = fact["subject"]
            fact_object = fact["object"]
            if not self.view.is_type(fact_subject, MOLECULAR_FUNCTION):
                continue
            if not self.view.is_type(fact_object, MOLECULAR_FUNCTION):
                continue

            subject_activities = self.activities_by_mf_id.get(fact_subject, [])
            object_activities = self.activities_by_mf_id.get(fact_object, [])
            if not subject_activities or not object_activities:
                continue
            if len(subject_activities) > 1:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.MULTIPLE_ACTIVITIES,
                        message=f"Multiple activities for subject {fact_subject}",
                        entity_id=fact_subject,
                    )
                )
            if len(object_activities) > 1:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.MULTIPLE_ACTIVITIES,
                        message=f"Multiple activities for object {fact_object}",
                        entity_id=fact_object,
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
        annotations = self.view.get_annotations(self.minerva_obj)
        annotations_mv = self.view.get_annotations_multivalued(self.minerva_obj)

        objects = [
            Object(
                id=obj["id"],
                label=obj.get("label"),
            )
            for obj in self.view.all_objects()
        ]

        provenance = ProvenanceInfo(
            contributor=annotations_mv.get("contributor"),
            date=annotations.get("date", None),
            provided_by=annotations_mv.get("providedBy"),
        )

        taxon, additional_taxa = self._resolve_model_taxa(annotations, annotations_mv)

        num_processed_facts = len(self.processed_facts)
        num_facts = len(self.view.all_facts())
        if num_processed_facts < num_facts:
            for fact in self.view.all_facts():
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

        if additional_taxa:
            cam.additional_taxa = additional_taxa

        return TranslationResult(result=cam, warnings=list(self.translation_warnings))

    def _resolve_model_taxa(
        self, annotations: dict, annotations_mv: dict
    ) -> tuple[str | None, list[str]]:
        """Resolve primary and additional taxa for the model."""
        all_taxa = annotations_mv.get(TaxonVocabulary.TAXON_ANNOTATION_KEY, [])
        legacy_taxon = annotations.get(TaxonVocabulary.LEGACY_TAXON_KEY)
        if legacy_taxon and legacy_taxon not in all_taxa:
            all_taxa.append(legacy_taxon)

        if not all_taxa:
            return None, []
        if len(all_taxa) == 1:
            return all_taxa[0], []

        host_matches = [t for t in all_taxa if TaxonVocabulary.is_host_taxon(t)]
        if host_matches:
            taxon = host_matches[0]
            additional_taxa = [t for t in all_taxa if t != taxon]
        else:
            taxon = all_taxa[0]
            additional_taxa = all_taxa[1:]

        return taxon, additional_taxa


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
