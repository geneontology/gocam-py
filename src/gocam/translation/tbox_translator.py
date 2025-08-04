from dataclasses import dataclass, field
from typing import Iterator, Union, Iterable, Optional, Dict, List

from mkdocs.commands.serve import serve
from oaklib.datamodels.vocabulary import IN_TAXON, RDFS_LABEL, PART_OF, IS_A, HAS_PART, MOLECULAR_FUNCTION, HAS_OUTPUT
from pyhornedowl.pyhornedowl import PyIndexedOntology

from gocam.datamodel import Model, Activity, TermAssociation, BiologicalProcessAssociation, \
    CellularAnatomicalEntityTermObject, CellularAnatomicalEntityAssociation, CellTypeAssociation, \
    GrossAnatomyAssociation, EnabledByProteinComplexAssociation, EnabledByGeneProductAssociation, MoleculeAssociation, \
    MolecularFunctionAssociation
from gocam.translation.minerva_wrapper import OCCURS_IN, HAS_INPUT, HAS_PRIMARY_INPUT, HAS_PRIMARY_OUTPUT, ENABLED_BY
import pyhornedowl
from pyhornedowl.model import SubClassOf, ObjectSomeValuesFrom, IRI, Class, ObjectProperty, DeclareClass, \
    AnnotationAssertion, Annotation, SimpleLiteral
import logging

logger = logging.getLogger(__name__)

HAPPENS_DURING = "RO:0002092"

SIMPLE_AXIOM = Union[SubClassOf, AnnotationAssertion]


@dataclass
class TBoxTranslator:
    """
    Translates Models to TBox OWL axioms.

    This translation is intended to capture faithfully OWL-level semantics of
    GO-CAM objects.
    """
    ontology: PyIndexedOntology = field(default_factory=lambda: PyIndexedOntology())
    _label_map: Optional[Dict[str, str]] = None

    def __post_init__(self):
        self.ontology = PyIndexedOntology()
        self.ontology.prefix_mapping.add_default_prefix_names()
        from gocam.datamodel.gocam import linkml_meta
        prefixes = linkml_meta.root["prefixes"]
        for pm in prefixes.values():
            p = pm.get("prefix_prefix")
            ref = pm.get("prefix_reference")
            self.ontology.add_prefix_mapping(p, ref)

    def _label(self, model: Model, id: str) -> str:
        if not self._label_map:
            self._label_map = {}
        if id not in self._label_map:
            # note: preserves map from last iteration...
            for term_obj in model.objects:
                self._label_map[term_obj.id] = term_obj.label
        return self._label_map.get(id, id)

    def load_models(self, models: Iterable[Model]) -> None:
        """
        Loads a collection of models into the current tbox ontology.

        :param models:
        :return:
        """
        for model in models:
            logger.info(f"Translating model {model.id} with title {model.title}")
            for axiom in self.translate_model_iter(model):
                self.ontology.add_axiom(axiom)

    def save_ontology(self, path: str, serialization : Optional[str] = None) -> None:
        """
        Save the current ontology to a file.

        :param path: The file path to save the ontology.
        :param serialization: The serialization format (e.g., "ofn", "turtle", "owlxml"). If None, defaults to "ofn".
        """
        if not serialization:
            serialization = "ofn"
        self.ontology.save_to_file(path, serialization=serialization)

    def translate_model_iter(self, model: Model) -> Iterator[SIMPLE_AXIOM]:
        """
        Translate an individual model into axioms

        :param model:
        :return:
        """
        yield from self.add_annotation(model.id, RDFS_LABEL, model.title)
        yield from self.add_edge(model.id, IS_A, "gocam:Model")
        for a in model.activities:
            yield from self.translate_activity(model, a)
            # yield from self.add_edge(a.id, PART_OF, model.id)
        taxon = model.taxon
        if taxon:
            yield from self.add_annotation(model.id, IN_TAXON, taxon)
        else:
            logger.warning(f"Model {model.id} has no taxon")

    def translate_activity(self, model: Model, activity: Activity) -> Iterator[SIMPLE_AXIOM]:
        if activity.molecular_function:
            mf_term_id = activity.molecular_function.term
            mf_label = self._label(model, mf_term_id)
        else:
            mf_term_id = MOLECULAR_FUNCTION
            mf_label = "UNKNOWN_FUNCTION"
        yield from self.add_edge(activity.id, IS_A, mf_term_id)
        yield from self.add_edge(activity.id, IS_A, "gocam:Activity")

        eb = activity.enabled_by
        if eb:
            if isinstance(eb, EnabledByProteinComplexAssociation):
                pc_id = self.make_id(model, eb.term, *eb.members)
                yield from self.add_edge(activity.id, ENABLED_BY, pc_id)
                member_labels = []
                for m in eb.members:
                    yield from self.add_edge(pc_id, HAS_PART, m)
                    member_labels.append(self._label(model, m))
                eb_label = f"{self._label(model, eb.term)}[{' '.join(member_labels)}]"
            elif isinstance(eb, EnabledByGeneProductAssociation):
                yield from self.add_edge(activity.id, ENABLED_BY, eb.term)
                eb_label = self._label(model, eb.term)
            else:
                raise NotImplementedError(f"Cannot handle {eb}")
        else:
            eb_label = "UNKNOWN"
        yield from self.add_annotation(activity.id, RDFS_LABEL, f"{eb_label} ({mf_label})")
        if activity.causal_associations:
            for ca in activity.causal_associations:
                predicate = ca.predicate
                yield from self.add_edge(activity.id, predicate, ca.downstream_activity)
        if activity.part_of:
            yield from self.translate_biological_process(model, activity.id, activity.part_of)
        else:
            yield from self.add_edge(activity.id, PART_OF, model.id)
        if activity.occurs_in:
            yield from self.translate_cellular_anatomical_entity(model, activity.id, activity.occurs_in)
        yield from self.translate_io(activity.id, activity.has_input, HAS_INPUT)
        yield from self.translate_io(activity.id, [activity.has_primary_input], HAS_PRIMARY_INPUT)
        yield from self.translate_io(activity.id, activity.has_output, HAS_OUTPUT)
        yield from self.translate_io(activity.id, [activity.has_primary_output], HAS_PRIMARY_OUTPUT)

    def _process_part_of_root(self, model: Model, ta: Union[Activity, BiologicalProcessAssociation]) -> str:
        if ta.part_of:
            return self._process_part_of_root(model, ta.part_of)
        elif isinstance(ta, Activity):
            return ta.id
        else:
            return ta.term


    def translate_biological_process(self, model: Model, child: str, ta: BiologicalProcessAssociation) -> Iterator[SIMPLE_AXIOM]:
        bp_term_id = ta.term
        bp_id = self.make_id(model, bp_term_id)
        label = f"{self._label(model, bp_term_id)} in {self._label(model, model.id)}"
        yield from self.add_annotation(bp_id, RDFS_LABEL, label)
        yield from self.add_edge(child, PART_OF, bp_id)
        yield from self.add_edge(bp_id, IS_A, bp_term_id)
        if ta.part_of:
            # recursive
            yield from self.translate_biological_process(model, bp_id, ta.part_of)
        else:
            yield from self.add_edge(bp_id, PART_OF, model.id)
        if ta.happens_during:
            yield from self.add_edge(bp_id, HAPPENS_DURING, ta.happens_during)

    def translate_cellular_anatomical_entity(self, model: Model, child: str, ta: CellularAnatomicalEntityAssociation,) -> Iterator[SIMPLE_AXIOM]:
        cc_term_id = ta.term
        cc_id = self.make_id(model, cc_term_id)
        label = f"{self._label(model, cc_term_id)} in {self._label(model, model.id)}"
        yield from self.add_annotation(cc_id, RDFS_LABEL, label)
        yield from self.add_edge(child, OCCURS_IN, cc_id)
        yield from self.add_edge(cc_id, IS_A, cc_term_id)
        if ta.part_of:
            # recursive; predicate switches to
            yield from self.translate_anatomical(model, cc_id, ta.part_of)

    def translate_anatomical(self, model: Model, child: str, ta: Union[CellTypeAssociation, GrossAnatomyAssociation]):
        anat_term_id = ta.term
        anat_id = self.make_id(model, anat_term_id)
        yield from self.add_edge(child, PART_OF, anat_id)
        yield from self.add_edge(anat_id, IS_A, anat_term_id)
        if ta.part_of:
            # recursive; predicate switches to
            yield from self.translate_anatomical(model, anat_id, ta.part_of)

    def translate_io(self, subject: str, tas: Optional[List[MoleculeAssociation]], predicate: str) -> Iterator[SubClassOf]:
        if tas:
            for ta in tas:
                if ta and ta.term and subject:
                    yield from self.add_edge(subject, predicate, ta.term)

    def add_annotation(self, subject: str, predicate: str, object: str) -> Iterator[AnnotationAssertion]:
        p_iri = self.iri(predicate)
        ap = self.ontology.annotation_property(str(p_iri))
        lv = SimpleLiteral(object)
        ann = Annotation(
            ap,
            lv,
        )
        yield AnnotationAssertion(self.iri(subject), ann)

    def add_edge(self, subject: str, predicate: str, object: str) -> Iterator[SubClassOf]:
        if predicate == IS_A:
            yield SubClassOf(self.get_class(subject), self.get_class(object))
        else:
            yield SubClassOf(self.get_class(subject),
                             ObjectSomeValuesFrom(self.get_predicate(predicate),
                                                  self.get_class(object)))

    def iri(self, id: str) -> IRI:
        if "/" in id:
            id = id.replace("/", "__")
        try:
            return self.ontology.curie(id)
        except ValueError as e:
            if ":" not in id:
                id = f"TEMP:{id}"
                logger.warning(f"Assuming {id} is a CURIE, but it does not contain a prefix")
            id_parts = id.split(":")
            pfx = id_parts[0]
            self.ontology.add_prefix_mapping(pfx, f"http://identifiers.org/{pfx}:")
            try:
                return self.ontology.curie(id)
            except ValueError as e:
                raise ValueError(f"{id} not a CURIE")

    def get_class(self, id: str) -> Class:
        iri = self.iri(id)
        return self.ontology.clazz(str(iri))

    def get_predicate(self, id: str) -> ObjectProperty:
        return self.ontology.object_property(str(self.iri(id)))

    def make_id(self, model: Model, *other_ids: str) -> str:
        other_ids_flat = "_".join(other_ids)
        return f"{model.id}_{other_ids_flat}"



