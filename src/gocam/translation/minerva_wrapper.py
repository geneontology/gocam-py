import logging
from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass, field

import requests
from pydantic import BaseModel, Field

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
    GrossAnatomyAssociation,
    Model,
    ModelStateEnum,
    MolecularFunctionAssociation,
    MoleculeAssociation,
    MoleculeNode,
    Object,
    PartOfProteinComplexAssociation,
    PhaseAssociation,
    ProteinComplexHasPartAssociation,
    ProvenanceInfo,
)
from gocam.translation.result import TranslationResult, TranslationWarning, WarningType
from gocam.vocabulary import (
    INFORMATION_BIOMACROMOLECULE,
    MOLECULAR_FUNCTION,
    PROTEIN_CONTAINING_COMPLEX,
    Relation,
    TaxonVocabulary,
)

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


class Annotation(BaseModel):
    key: str
    value: str
    value_type: str | None = Field(default=None, alias="value-type")


class IndividualType(BaseModel):
    type: str
    id: str = ""
    label: str | None = None


class Annotated(BaseModel):
    annotations: list[Annotation] = Field(default_factory=list)


class Individual(Annotated):
    id: str
    type: list[IndividualType] = Field(default_factory=list)
    root_type: list[IndividualType] = Field(default_factory=list, alias="root-type")


class Fact(Annotated):
    subject: str
    property: str
    object: str


class MinervaObject(Annotated):
    id: str
    individuals: list[Individual] = Field(default_factory=list)
    facts: list[Fact] = Field(default_factory=list)


@dataclass
class MinervaView:
    """
    Indexed, read-only view of raw Minerva JSON data.
    """

    raw_json: MinervaObject
    _facts_by_property: defaultdict[str, list[Fact]] = field(
        init=False, default_factory=lambda: defaultdict(list)
    )
    _facts_by_subject_property: defaultdict[tuple[str, str], list[Fact]] = field(
        init=False, default_factory=lambda: defaultdict(list)
    )
    _facts_by_object_property: defaultdict[tuple[str, str], list[Fact]] = field(
        init=False, default_factory=lambda: defaultdict(list)
    )
    _individual_id_to_term: dict[str, str] = field(init=False, default_factory=dict)
    _individual_id_to_root_types: dict[str, list[IndividualType]] = field(
        init=False, default_factory=dict
    )
    _individual_id_to_annotations: dict[str, dict[str, str]] = field(
        init=False, default_factory=dict
    )
    _individual_id_to_annotations_multivalued: dict[str, dict[str, list[str]]] = field(
        init=False, default_factory=dict
    )
    _objects_by_id: dict[str, IndividualType] = field(init=False, default_factory=dict)

    def __post_init__(self):
        self._facts_by_property = defaultdict(list)
        self._facts_by_subject_property = defaultdict(list)
        self._facts_by_object_property = defaultdict(list)

        for individual in self.raw_json.individuals:
            ind_id = individual.id
            self._individual_id_to_root_types[ind_id] = [
                x for x in individual.root_type if x
            ]
            self._individual_id_to_annotations[ind_id] = self._extract_annotations(
                individual
            )
            self._individual_id_to_annotations_multivalued[ind_id] = (
                self._extract_annotations_multivalued(individual)
            )

            for type_ in individual.type:
                if type_.type == "complement":
                    continue
                type_id = type_.id
                if type_id:
                    self._objects_by_id[type_id] = type_
                    self._individual_id_to_term[ind_id] = type_id

        for fact in self.raw_json.facts:
            self._facts_by_property[fact.property].append(fact)
            self._facts_by_subject_property[(fact.subject, fact.property)].append(fact)
            self._facts_by_object_property[(fact.object, fact.property)].append(fact)

    def _normalize_property(self, prop: str) -> str:
        """Normalize a property URI."""
        if "/" in prop:
            return prop.split("/")[-1]
        return prop

    def _extract_annotations(self, obj: Annotated) -> dict[str, str]:
        """Extract single-valued annotations."""
        return {self._normalize_property(a.key): a.value for a in obj.annotations}

    def _extract_annotations_multivalued(self, obj: Annotated) -> dict[str, list[str]]:
        """Extract multi-valued annotations."""
        anns = defaultdict(list)
        for a in obj.annotations:
            key = self._normalize_property(a.key)
            value = a.value
            anns[key].append(value)
        return anns

    def get_term(self, individual_id: str) -> str | None:
        """Get the term (GO, CHEBI, ECO, etc.) associated with an individual."""
        return self._individual_id_to_term.get(individual_id)

    def get_root_types(self, individual_id: str) -> list[str]:
        """Get the root types associated with an individual."""
        return [
            rt.id for rt in self._individual_id_to_root_types.get(individual_id, [])
        ]

    def is_type(self, individual_id: str, type_uri: str) -> bool:
        """Check if an individual has a specific root type."""
        return type_uri in self.get_root_types(individual_id)

    def get_annotations(self, id_or_dict: str | Annotated) -> dict[str, str]:
        """Get single-valued annotations for an individual (by ID) or a fact (by dict)."""
        if isinstance(id_or_dict, str):
            return self._individual_id_to_annotations.get(id_or_dict, {})
        return self._extract_annotations(id_or_dict)

    def get_annotations_multivalued(
        self, id_or_dict: str | Annotated
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
    ) -> list[Fact]:
        """Query facts by subject, object, and/or predicate."""
        if subject and property and not object:
            return self._facts_by_subject_property.get((subject, property), [])
        if object and property and not subject:
            return self._facts_by_object_property.get((object, property), [])
        if property and not subject and not object:
            return self._facts_by_property.get(property, [])
        # Fallback to iteration for other filters
        return [
            f
            for f in self.raw_json.facts
            if (not subject or f.subject == subject)
            and (not object or f.object == object)
            and (not property or f.property == property)
        ]

    def get_individual_label(self, individual_id: str) -> str:
        """Get the string representation of an individual, if it exists."""
        term_id = self._individual_id_to_term.get(individual_id)
        if term_id is None:
            return individual_id
        obj = self._objects_by_id.get(term_id)
        if obj is None:
            return term_id
        label = obj.label
        formatted = f"{label} [{term_id}]" if label else term_id
        types = self._individual_id_to_root_types.get(individual_id, [])
        if types and types[0].label:
            # Show the first one because it is the most relevant one
            formatted += f" ({types[0].label})"
        return formatted

    def get_prop_label(self, prop: str):
        """Get the string representation of a property, if it exists."""
        try:
            relation = Relation(prop)
            return relation.name.lower().replace("_", " ")
        except ValueError:
            return prop

    def all_objects(self) -> list[IndividualType]:
        """Get all object definitions."""
        return list(self._objects_by_id.values())

    def all_facts(self) -> list[Fact]:
        """Get all facts."""
        return self.raw_json.facts


@dataclass
class MinervaTranslator:
    """
    Translates a Minerva JSON object to a GO-CAM Model.
    """

    minerva_obj: MinervaObject
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

    def _add_skipped_fact_warning(self, fact: Fact):
        """Add a warning that a fact was skipped because the relevant attribute on the Activity already had a value."""
        prop_name = Relation.get_name(fact.property)
        self.translation_warnings.add(
            TranslationWarning(
                type=WarningType.SKIPPED_FACT,
                message=f"'{self.view.get_individual_label(fact.subject)}' already has '{prop_name}' relationship; "
                f"ignoring '{prop_name}' relationship to '{self.view.get_individual_label(fact.object)}'",
                entity_id=fact.subject,
            )
        )

    def _add_missing_term_warning(self, individual_id: str, fact: Fact):
        """Add a warning that an individual is missing a term, including context about where the issue was encountered."""
        bad_individual = "subject" if fact.subject == individual_id else "object"
        self.translation_warnings.add(
            TranslationWarning(
                type=WarningType.MISSING_TERM,
                message=f"While processing '{self.view.get_individual_label(fact.subject)}' "
                f"{Relation.get_name(fact.property)} '{self.view.get_individual_label(fact.object)}', "
                f"the {bad_individual} could not be mapped to a term; skipping this fact.",
                entity_id=individual_id,
            )
        )

    def _process_fact(self, fact: Fact) -> tuple[list[EvidenceItem], ProvenanceInfo]:
        """Process a fact to extract evidence and provenance information.

        This method also marks the fact as processed by adding it to the processed_facts set, which
        is used to track which facts have been handled during translation and identify any that
        were not.
        """
        self.processed_facts.add((fact.subject, fact.property, fact.object))
        return self._extract_evidence_and_provenance(fact)

    def _extract_evidence_and_provenance(
        self, id_or_dict: str | Fact
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
                    entity_id=self.minerva_obj.id,
                )
            )
        for enabled_by_fact in enabled_by_facts:
            evs, prov = self._process_fact(enabled_by_fact)
            enabled_by_subject = enabled_by_fact.subject
            enabled_by_object = enabled_by_fact.object
            subject_term = self.view.get_term(enabled_by_subject)
            object_term = self.view.get_term(enabled_by_object)

            if subject_term is None:
                self._add_missing_term_warning(enabled_by_subject, enabled_by_fact)
                continue
            if object_term is None:
                self._add_missing_term_warning(enabled_by_object, enabled_by_fact)
                continue

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
            individual_label = self.view.get_individual_label(individual_id)
            self.translation_warnings.add(
                TranslationWarning(
                    type=WarningType.UNKNOWN_ENABLED_BY_TYPE,
                    message=f"'{individual_label}' is used in enabled by relationship, but it is not "
                    f"of type protein-containing complex [GO:0032991] or information biomacromolecule "
                    f"[CHEBI:33695]; assuming information biomacromolecule / gene product.",
                    entity_id=individual_id,
                )
            )

        return self._build_gene_product_association(individual_id, term, evidence, prov)

    def _build_protein_complex_association(
        self,
        individual_id: str,
        term: str,
        evidence: list[EvidenceItem],
        prov: ProvenanceInfo,
    ) -> EnabledByProteinComplexAssociation:
        """Build an EnabledByProteinComplexAssociation by gathering member facts."""

        enabled_by_association = EnabledByProteinComplexAssociation(
            term=term,
            evidence=evidence,
            provenances=[prov],
        )

        for has_part_fact in self.view.get_facts(
            subject=individual_id, property=Relation.HAS_PART
        ):
            has_part_association = self._build_protein_complex_has_part_association(
                has_part_fact
            )
            if has_part_association is None:
                continue
            if enabled_by_association.has_part is None:
                enabled_by_association.has_part = []
            enabled_by_association.has_part.append(has_part_association)

        return enabled_by_association

    def _build_gene_product_association(
        self,
        individual_id: str,
        term: str,
        evidence: list[EvidenceItem],
        prov: ProvenanceInfo,
    ) -> EnabledByGeneProductAssociation:
        """Build an EnabledByGeneProductAssociation."""

        enabled_by_association = EnabledByGeneProductAssociation(
            term=term,
            evidence=evidence,
            provenances=[prov],
        )

        for part_of_fact in self.view.get_facts(
            subject=individual_id, property=Relation.PART_OF
        ):
            part_of_association = self._build_part_of_protein_complex_association(
                part_of_fact
            )
            if part_of_association is None:
                continue
            if enabled_by_association.part_of is None:
                enabled_by_association.part_of = []
            enabled_by_association.part_of.append(part_of_association)

        return enabled_by_association

    def _build_protein_complex_has_part_association(
        self, fact: Fact, visited: frozenset[str] | None = None
    ) -> ProteinComplexHasPartAssociation | None:
        """Build a ProteinComplexHasPartAssociation from a has_part fact."""
        if visited is None:
            visited = frozenset()

        evs, prov = self._process_fact(fact)
        part_id = fact.object
        part_term = self.view.get_term(part_id)

        if part_term is None:
            self._add_missing_term_warning(part_id, fact)
            return None

        has_part_association = ProteinComplexHasPartAssociation(
            term=part_term,
            evidence=evs,
            provenances=[prov],
        )

        for part_of_fact in self.view.get_facts(
            subject=part_id, property=Relation.PART_OF
        ):
            part_of_object = part_of_fact.object
            if part_of_object in visited:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.INVALID_CIRCULAR_RELATIONSHIP,
                        message=f"Circular relationship detected for {self.view.get_individual_label(part_of_object)}",
                        entity_id=part_of_object,
                    )
                )
                continue
            part_of_association = self._build_part_of_protein_complex_association(
                part_of_fact, visited=visited.union({part_of_object})
            )
            if part_of_association is None:
                continue
            if has_part_association.part_of is None:
                has_part_association.part_of = []
            has_part_association.part_of.append(part_of_association)

        return has_part_association

    def _build_part_of_protein_complex_association(
        self, fact: Fact, visited: frozenset[str] | None = None
    ) -> PartOfProteinComplexAssociation | None:
        """Build a PartOfProteinComplexAssociation from a part_of fact."""
        if visited is None:
            visited = frozenset()

        evs, prov = self._process_fact(fact)
        complex_id = fact.object
        complex_term = self.view.get_term(complex_id)
        if complex_term is None:
            self._add_missing_term_warning(complex_id, fact)
            return None

        part_of_association = PartOfProteinComplexAssociation(
            term=complex_term,
            evidence=evs,
            provenances=[prov],
        )

        for has_part_fact in self.view.get_facts(
            subject=complex_id, property=Relation.HAS_PART
        ):
            has_part_object = has_part_fact.object
            if has_part_object in visited:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.INVALID_CIRCULAR_RELATIONSHIP,
                        message=f"Circular relationship detected for {self.view.get_individual_label(has_part_object)}",
                        entity_id=has_part_object,
                    )
                )
                continue
            has_part_association = self._build_protein_complex_has_part_association(
                has_part_fact, visited=visited.union({has_part_object})
            )
            if has_part_association is None:
                continue
            if part_of_association.has_part is None:
                part_of_association.has_part = []
            part_of_association.has_part.append(has_part_association)

        return part_of_association

    def _add_molecule_node(self, individual_id: str, fact: Fact) -> MoleculeNode | None:
        """Add a molecule node for a given individual ID, if it doesn't already exist."""
        if individual_id in self.molecule_nodes_by_id:
            return self.molecule_nodes_by_id[individual_id]

        term = self.view.get_term(individual_id)
        if term is None:
            self._add_missing_term_warning(individual_id, fact)
            return None

        molecule_node = MoleculeNode(
            id=individual_id,
            term=term,
        )

        for located_in_fact in self.view.get_facts(
            subject=individual_id, property=Relation.LOCATED_IN
        ):
            association = self._build_cellular_anatomical_entity_association(
                located_in_fact
            )
            if association is None:
                continue
            if molecule_node.located_in is None:
                molecule_node.located_in = association
            else:
                self._add_skipped_fact_warning(located_in_fact)

        self.molecule_nodes_by_id[individual_id] = molecule_node
        return molecule_node

    def _add_molecule_association(
        self, fact: Fact, *, is_inverse: bool = False
    ) -> None:
        """Add a molecule association to the relevant activity based on a fact."""
        if not is_inverse:
            mf_individual_id = fact.subject
            molecule_individual_id = fact.object
            predicate = fact.property
        else:
            mf_individual_id = fact.object
            molecule_individual_id = fact.subject
            try:
                predicate_relation = Relation(fact.property)
                predicate = MOLECULAR_ASSOCIATION_INVERSE_PROPERTIES[predicate_relation]
            except (ValueError, KeyError):
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.UNKNOWN_PROPERTY_INVERSE,
                        message=f"Fact property '{fact.property}' does not have a known inverse",
                        entity_id=fact.property,
                    )
                )
                return

        molecule_node = self._add_molecule_node(molecule_individual_id, fact)
        if molecule_node is None:
            return

        activities = self.activities_by_mf_id.get(mf_individual_id)
        if not activities:
            subj_str = self.view.get_individual_label(fact.subject)
            obj_str = self.view.get_individual_label(fact.object)
            prop_str = Relation.get_name(fact.property)
            self.translation_warnings.add(
                TranslationWarning(
                    type=WarningType.SKIPPED_FACT,
                    message=f"Subject '{subj_str}' has '{prop_str}' relationship to '{obj_str}', but that subject is not the molecular function of an activity unit; skipping molecular association",
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
        self, fact: Fact, visited: frozenset[str] | None = None
    ) -> BiologicalProcessAssociation | None:
        """Recursively build a BiologicalProcessAssociation."""
        evidence, provenance = self._process_fact(fact)
        individual_id = fact.object

        term = self.view.get_term(individual_id)
        if term is None:
            self._add_missing_term_warning(individual_id, fact)
            return None

        if visited is None:
            visited = frozenset()
        if individual_id in visited:
            self.translation_warnings.add(
                TranslationWarning(
                    type=WarningType.INVALID_CIRCULAR_RELATIONSHIP,
                    message=f"Circular relationship detected for {self.view.get_individual_label(individual_id)}",
                    entity_id=individual_id,
                )
            )
            return None
        next_visited = visited.union({individual_id})

        association = BiologicalProcessAssociation(
            term=term,
            evidence=evidence,
            provenances=[provenance],
        )

        for happens_during_fact in self.view.get_facts(
            subject=individual_id, property=Relation.HAPPENS_DURING
        ):
            phase_association = self._build_phase_association(happens_during_fact)
            if phase_association is None:
                continue
            if association.happens_during is None:
                association.happens_during = phase_association
            else:
                self._add_skipped_fact_warning(happens_during_fact)

        for part_of_fact in self.view.get_facts(
            subject=individual_id, property=Relation.PART_OF
        ):
            bp_association = self._build_biological_process_association(
                part_of_fact, visited=next_visited
            )
            if bp_association is None:
                continue
            if association.part_of is None:
                association.part_of = bp_association
            else:
                self._add_skipped_fact_warning(part_of_fact)

        return association

    def _build_phase_association(self, fact: Fact) -> PhaseAssociation | None:
        """Build a PhaseAssociation."""
        evidence, provenance = self._process_fact(fact)
        individual_id = fact.object
        term = self.view.get_term(individual_id)
        if term is None:
            self._add_missing_term_warning(individual_id, fact)
            return None
        return PhaseAssociation(
            term=term,
            evidence=evidence,
            provenances=[provenance],
        )

    def _build_cellular_anatomical_entity_association(
        self, fact: Fact
    ) -> CellularAnatomicalEntityAssociation | None:
        """Build a CellularAnatomicalEntityAssociation."""
        evidence, provenance = self._process_fact(fact)
        individual_id = fact.object
        term = self.view.get_term(individual_id)
        if term is None:
            self._add_missing_term_warning(individual_id, fact)
            return None
        association = CellularAnatomicalEntityAssociation(
            term=term,
            evidence=evidence,
            provenances=[provenance],
        )

        for part_of_fact in self.view.get_facts(
            subject=individual_id, property=Relation.PART_OF
        ):
            cell_type_association = self._build_cell_type_association(part_of_fact)
            if cell_type_association is None:
                continue
            if association.part_of is None:
                association.part_of = cell_type_association
            else:
                self._add_skipped_fact_warning(part_of_fact)

        return association

    def _build_cell_type_association(self, fact: Fact) -> CellTypeAssociation | None:
        """Build a CellTypeAssociation."""
        evidence, provenance = self._process_fact(fact)
        individual_id = fact.object
        term = self.view.get_term(individual_id)
        if term is None:
            self._add_missing_term_warning(individual_id, fact)
            return None
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
            if gross_anatomy_association is None:
                continue
            if association.part_of is None:
                association.part_of = gross_anatomy_association
            else:
                self._add_skipped_fact_warning(part_of_fact)

        return association

    def _build_gross_anatomy_association(
        self, fact: Fact, visited: frozenset[str] | None = None
    ) -> GrossAnatomyAssociation | None:
        """Recursively build a GrossAnatomyAssociation."""
        evidence, provenance = self._process_fact(fact)
        individual_id = fact.object

        if visited is None:
            visited = frozenset()
        if individual_id in visited:
            self.translation_warnings.add(
                TranslationWarning(
                    type=WarningType.INVALID_CIRCULAR_RELATIONSHIP,
                    message=f"Circular relationship detected for {self.view.get_individual_label(individual_id)}",
                    entity_id=individual_id,
                )
            )
            return None
        next_visited = visited.union({individual_id})

        term = self.view.get_term(individual_id)
        if term is None:
            self._add_missing_term_warning(individual_id, fact)
            return None

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
            if gross_anatomy_association is None:
                continue
            if association.part_of is None:
                association.part_of = gross_anatomy_association
            else:
                self._add_skipped_fact_warning(part_of_fact)

        return association

    def _process_biological_process_associations(self):
        """Process part_of facts to build biological process associations."""
        for fact in self.view.get_facts(property=Relation.PART_OF):
            part_of_subject = fact.subject

            activities = self.activities_by_mf_id.get(part_of_subject)
            if not activities:
                continue

            association = self._build_biological_process_association(fact)
            if association is None:
                continue

            for activity in activities:
                if activity.part_of is None:
                    activity.part_of = association
                else:
                    self._add_skipped_fact_warning(fact)

    def _process_occurs_in_associations(self):
        """Process occurs_in facts to build cellular anatomical entity associations."""
        for fact in self.view.get_facts(property=Relation.OCCURS_IN):
            occurs_in_subject = fact.subject

            activities = self.activities_by_mf_id.get(occurs_in_subject)
            if not activities:
                continue

            association = self._build_cellular_anatomical_entity_association(fact)
            if association is None:
                continue

            for activity in activities:
                if activity.occurs_in is None:
                    activity.occurs_in = association
                else:
                    self._add_skipped_fact_warning(fact)

    def _process_happens_during_associations(self):
        """Process happens_during facts to build phase associations."""
        for fact in self.view.get_facts(property=Relation.HAPPENS_DURING):
            happens_during_subject = fact.subject

            activities = self.activities_by_mf_id.get(happens_during_subject)
            if not activities:
                continue

            association = self._build_phase_association(fact)
            if association is None:
                continue

            for activity in activities:
                if activity.happens_during is None:
                    activity.happens_during = association
                else:
                    self._add_skipped_fact_warning(fact)

    def _process_molecule_associations(self):
        """Process facts representing associations between molecular functions and molecules."""
        for fact in self.view.all_facts():
            fact_property = fact.property
            if fact_property in MOLECULAR_ASSOCIATION_PROPERTIES:
                self._add_molecule_association(fact)
            elif fact_property in MOLECULAR_ASSOCIATION_INVERSE_PROPERTIES:
                self._add_molecule_association(fact, is_inverse=True)

    def _process_causal_associations(self):
        """Process facts representing causal associations between activities."""
        for fact in self.view.all_facts():
            fact_property = fact.property
            fact_subject = fact.subject
            fact_object = fact.object
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
                        message=f"Multiple activities for subject {self.view.get_individual_label(fact_subject)}",
                        entity_id=fact_subject,
                    )
                )
            if len(object_activities) > 1:
                self.translation_warnings.add(
                    TranslationWarning(
                        type=WarningType.MULTIPLE_ACTIVITIES,
                        message=f"Multiple activities for object {self.view.get_individual_label(fact_object)}",
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
                id=obj.id,
                label=obj.label,
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
                fact_tuple = (fact.subject, fact.property, fact.object)
                if fact_tuple not in self.processed_facts:
                    subj = self.view.get_individual_label(fact.subject)
                    obj = self.view.get_individual_label(fact.object)
                    prop = Relation.get_name(fact.property)
                    self.translation_warnings.add(
                        TranslationWarning(
                            type=WarningType.UNHANDLED_FACT,
                            message=f"Fact not representable in GO-CAM model: {subj} —{prop}→ {obj}",
                            entity_id=self.minerva_obj.id,
                        )
                    )

        cam = Model(
            id=self.minerva_obj.id,
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
                        message=f"Invalid model state: {state_annotation}",
                        entity_id=self.minerva_obj.id,
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
        minerva_model = MinervaObject.model_validate(minerva_obj)
        translator = MinervaTranslator(minerva_model)
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
