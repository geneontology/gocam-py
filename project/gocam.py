# Auto generated from gocam.yaml by pythongen.py version: 0.0.1
# Generation date: 2025-07-29T09:39:31
# Schema: gocam
#
# id: https://w3id.org/gocam
# description: Gene Ontology Causal Activity Model (GO-CAM) Schema.
#
#   This schema provides a way of representing causal pathway [Models](Model.md). A model consists of a set of
#   [Activity](Activity.md) objects, where each activity object represents the function of either an [individual
#   gene product](EnabledByGeneProductAssociation), a [protein complex of gene products](EnabledByGeneProductAssociation),
#   or a set of possible gene products.
#
#   Each [Models](Model.md) has associated metadata slots. Some slots such as [id](id.md), [title](title.md),
#   and [status](status.md) are *required*.
# license: https://creativecommons.org/publicdomain/zero/1.0/

import dataclasses
import re
from dataclasses import dataclass
from datetime import (
    date,
    datetime,
    time
)
from typing import (
    Any,
    ClassVar,
    Dict,
    List,
    Optional,
    Union
)

from jsonasobj2 import (
    JsonObj,
    as_dict
)
from linkml_runtime.linkml_model.meta import (
    EnumDefinition,
    PermissibleValue,
    PvFormulaOptions
)
from linkml_runtime.utils.curienamespace import CurieNamespace
from linkml_runtime.utils.enumerations import EnumDefinitionImpl
from linkml_runtime.utils.formatutils import (
    camelcase,
    sfx,
    underscore
)
from linkml_runtime.utils.metamodelcore import (
    bnode,
    empty_dict,
    empty_list
)
from linkml_runtime.utils.slot import Slot
from linkml_runtime.utils.yamlutils import (
    YAMLRoot,
    extended_float,
    extended_int,
    extended_str
)
from rdflib import (
    Namespace,
    URIRef
)

from linkml_runtime.linkml_model.types import Boolean, Integer, String, Uriorcurie
from linkml_runtime.utils.metamodelcore import Bool, URIorCURIE

metamodel_version = "1.7.0"
version = None

# Namespaces
BFO = CurieNamespace('BFO', 'http://purl.obolibrary.org/obo/BFO_')
CHEBI = CurieNamespace('CHEBI', 'http://example.org/UNKNOWN/CHEBI/')
CL = CurieNamespace('CL', 'http://example.org/UNKNOWN/CL/')
DDANAT = CurieNamespace('DDANAT', 'http://example.org/UNKNOWN/DDANAT/')
DOI = CurieNamespace('DOI', 'http://example.org/UNKNOWN/DOI/')
ECO = CurieNamespace('ECO', 'http://purl.obolibrary.org/obo/ECO_')
FAO = CurieNamespace('FAO', 'http://example.org/UNKNOWN/FAO/')
GO = CurieNamespace('GO', 'http://purl.obolibrary.org/obo/GO_')
GOREF = CurieNamespace('GOREF', 'http://example.org/UNKNOWN/GOREF/')
NCBITAXON = CurieNamespace('NCBITaxon', 'http://purl.obolibrary.org/obo/NCBITaxon_')
OBAN = CurieNamespace('OBAN', 'http://purl.org/oban/')
PMID = CurieNamespace('PMID', 'http://identifiers.org/pubmed/')
PO = CurieNamespace('PO', 'http://example.org/UNKNOWN/PO/')
RHEA = CurieNamespace('RHEA', 'http://rdf.rhea-db.org/')
RO = CurieNamespace('RO', 'http://purl.obolibrary.org/obo/RO_')
UBERON = CurieNamespace('UBERON', 'http://example.org/UNKNOWN/UBERON/')
UNIPROTKB = CurieNamespace('UniProtKB', 'http://purl.uniprot.org/uniprot/')
BIOLINK = CurieNamespace('biolink', 'https://w3id.org/biolink/vocab/')
DCE = CurieNamespace('dce', 'http://purl.org/dc/elements/1.1/')
DCT = CurieNamespace('dct', 'http://example.org/UNKNOWN/dct/')
DCTERMS = CurieNamespace('dcterms', 'http://purl.org/dc/terms/')
GOCAM = CurieNamespace('gocam', 'https://w3id.org/gocam/')
GOMODEL = CurieNamespace('gomodel', 'http://model.geneontology.org/')
GOSHAPES = CurieNamespace('goshapes', 'http://purl.obolibrary.org/obo/go/shapes/')
LEGO = CurieNamespace('lego', 'http://geneontology.org/lego/')
LINKML = CurieNamespace('linkml', 'https://w3id.org/linkml/')
OIO = CurieNamespace('oio', 'http://www.geneontology.org/formats/oboInOwl#')
ORCID = CurieNamespace('orcid', 'https://orcid.org/')
PAV = CurieNamespace('pav', 'http://purl.org/pav/')
RDFS = CurieNamespace('rdfs', 'http://example.org/UNKNOWN/rdfs/')
DEFAULT_ = GOCAM


# Types

# Class references
class ModelId(URIorCURIE):
    pass


class ActivityId(URIorCURIE):
    pass


class ObjectId(URIorCURIE):
    pass


class TermObjectId(ObjectId):
    pass


class PublicationObjectId(ObjectId):
    pass


class EvidenceTermObjectId(TermObjectId):
    pass


class MolecularFunctionTermObjectId(TermObjectId):
    pass


class BiologicalProcessTermObjectId(TermObjectId):
    pass


class CellularAnatomicalEntityTermObjectId(TermObjectId):
    pass


class MoleculeTermObjectId(TermObjectId):
    pass


class CellTypeTermObjectId(TermObjectId):
    pass


class GrossAnatomicalStructureTermObjectId(TermObjectId):
    pass


class PhaseTermObjectId(TermObjectId):
    pass


class InformationBiomacromoleculeTermObjectId(TermObjectId):
    pass


class GeneProductTermObjectId(InformationBiomacromoleculeTermObjectId):
    pass


class ProteinComplexTermObjectId(InformationBiomacromoleculeTermObjectId):
    pass


class TaxonTermObjectId(TermObjectId):
    pass


class PredicateTermObjectId(TermObjectId):
    pass


@dataclass(repr=False)
class Model(YAMLRoot):
    """
    A model of a biological program consisting of a set of causally connected activities.
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["Model"]
    class_class_curie: ClassVar[str] = "gocam:Model"
    class_name: ClassVar[str] = "Model"
    class_model_uri: ClassVar[URIRef] = GOCAM.Model

    id: Union[str, ModelId] = None
    title: str = None
    taxon: Optional[Union[str, TaxonTermObjectId]] = None
    additional_taxa: Optional[Union[Union[str, TaxonTermObjectId], list[Union[str, TaxonTermObjectId]]]] = empty_list()
    status: Optional[Union[str, "ModelStateEnum"]] = None
    comments: Optional[Union[str, list[str]]] = empty_list()
    activities: Optional[Union[dict[Union[str, ActivityId], Union[dict, "Activity"]], list[Union[dict, "Activity"]]]] = empty_dict()
    objects: Optional[Union[dict[Union[str, ObjectId], Union[dict, "Object"]], list[Union[dict, "Object"]]]] = empty_dict()
    provenances: Optional[Union[Union[dict, "ProvenanceInfo"], list[Union[dict, "ProvenanceInfo"]]]] = empty_list()
    query_index: Optional[Union[dict, "QueryIndex"]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, ModelId):
            self.id = ModelId(self.id)

        if self._is_empty(self.title):
            self.MissingRequiredField("title")
        if not isinstance(self.title, str):
            self.title = str(self.title)

        if self.taxon is not None and not isinstance(self.taxon, TaxonTermObjectId):
            self.taxon = TaxonTermObjectId(self.taxon)

        if not isinstance(self.additional_taxa, list):
            self.additional_taxa = [self.additional_taxa] if self.additional_taxa is not None else []
        self.additional_taxa = [v if isinstance(v, TaxonTermObjectId) else TaxonTermObjectId(v) for v in self.additional_taxa]

        if self.status is not None and not isinstance(self.status, ModelStateEnum):
            self.status = ModelStateEnum(self.status)

        if not isinstance(self.comments, list):
            self.comments = [self.comments] if self.comments is not None else []
        self.comments = [v if isinstance(v, str) else str(v) for v in self.comments]

        self._normalize_inlined_as_list(slot_name="activities", slot_type=Activity, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="objects", slot_type=Object, key_name="id", keyed=True)

        if not isinstance(self.provenances, list):
            self.provenances = [self.provenances] if self.provenances is not None else []
        self.provenances = [v if isinstance(v, ProvenanceInfo) else ProvenanceInfo(**as_dict(v)) for v in self.provenances]

        if self.query_index is not None and not isinstance(self.query_index, QueryIndex):
            self.query_index = QueryIndex(**as_dict(self.query_index))

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class Activity(YAMLRoot):
    """
    An individual activity in a causal model, representing the individual molecular activity of a single gene product
    or complex in the context of a particular model
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["Activity"]
    class_class_curie: ClassVar[str] = "gocam:Activity"
    class_name: ClassVar[str] = "Activity"
    class_model_uri: ClassVar[URIRef] = GOCAM.Activity

    id: Union[str, ActivityId] = None
    enabled_by: Optional[Union[dict, "EnabledByAssociation"]] = None
    molecular_function: Optional[Union[dict, "MolecularFunctionAssociation"]] = None
    occurs_in: Optional[Union[dict, "CellularAnatomicalEntityAssociation"]] = None
    part_of: Optional[Union[dict, "BiologicalProcessAssociation"]] = None
    has_input: Optional[Union[Union[dict, "MoleculeAssociation"], list[Union[dict, "MoleculeAssociation"]]]] = empty_list()
    has_primary_input: Optional[Union[dict, "MoleculeAssociation"]] = None
    has_output: Optional[Union[Union[dict, "MoleculeAssociation"], list[Union[dict, "MoleculeAssociation"]]]] = empty_list()
    has_primary_output: Optional[Union[dict, "MoleculeAssociation"]] = None
    causal_associations: Optional[Union[Union[dict, "CausalAssociation"], list[Union[dict, "CausalAssociation"]]]] = empty_list()
    provenances: Optional[Union[Union[dict, "ProvenanceInfo"], list[Union[dict, "ProvenanceInfo"]]]] = empty_list()

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, ActivityId):
            self.id = ActivityId(self.id)

        if self.enabled_by is not None and not isinstance(self.enabled_by, EnabledByAssociation):
            self.enabled_by = EnabledByAssociation(**as_dict(self.enabled_by))

        if self.molecular_function is not None and not isinstance(self.molecular_function, MolecularFunctionAssociation):
            self.molecular_function = MolecularFunctionAssociation(**as_dict(self.molecular_function))

        if self.occurs_in is not None and not isinstance(self.occurs_in, CellularAnatomicalEntityAssociation):
            self.occurs_in = CellularAnatomicalEntityAssociation(**as_dict(self.occurs_in))

        if self.part_of is not None and not isinstance(self.part_of, BiologicalProcessAssociation):
            self.part_of = BiologicalProcessAssociation(**as_dict(self.part_of))

        if not isinstance(self.has_input, list):
            self.has_input = [self.has_input] if self.has_input is not None else []
        self.has_input = [v if isinstance(v, MoleculeAssociation) else MoleculeAssociation(**as_dict(v)) for v in self.has_input]

        if self.has_primary_input is not None and not isinstance(self.has_primary_input, MoleculeAssociation):
            self.has_primary_input = MoleculeAssociation(**as_dict(self.has_primary_input))

        if not isinstance(self.has_output, list):
            self.has_output = [self.has_output] if self.has_output is not None else []
        self.has_output = [v if isinstance(v, MoleculeAssociation) else MoleculeAssociation(**as_dict(v)) for v in self.has_output]

        if self.has_primary_output is not None and not isinstance(self.has_primary_output, MoleculeAssociation):
            self.has_primary_output = MoleculeAssociation(**as_dict(self.has_primary_output))

        if not isinstance(self.causal_associations, list):
            self.causal_associations = [self.causal_associations] if self.causal_associations is not None else []
        self.causal_associations = [v if isinstance(v, CausalAssociation) else CausalAssociation(**as_dict(v)) for v in self.causal_associations]

        if not isinstance(self.provenances, list):
            self.provenances = [self.provenances] if self.provenances is not None else []
        self.provenances = [v if isinstance(v, ProvenanceInfo) else ProvenanceInfo(**as_dict(v)) for v in self.provenances]

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class EvidenceItem(YAMLRoot):
    """
    An individual piece of evidence that is associated with an assertion in a model
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["EvidenceItem"]
    class_class_curie: ClassVar[str] = "gocam:EvidenceItem"
    class_name: ClassVar[str] = "EvidenceItem"
    class_model_uri: ClassVar[URIRef] = GOCAM.EvidenceItem

    term: Optional[Union[str, EvidenceTermObjectId]] = None
    reference: Optional[Union[str, PublicationObjectId]] = None
    with_objects: Optional[Union[Union[str, ObjectId], list[Union[str, ObjectId]]]] = empty_list()
    provenances: Optional[Union[Union[dict, "ProvenanceInfo"], list[Union[dict, "ProvenanceInfo"]]]] = empty_list()

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.term is not None and not isinstance(self.term, EvidenceTermObjectId):
            self.term = EvidenceTermObjectId(self.term)

        if self.reference is not None and not isinstance(self.reference, PublicationObjectId):
            self.reference = PublicationObjectId(self.reference)

        if not isinstance(self.with_objects, list):
            self.with_objects = [self.with_objects] if self.with_objects is not None else []
        self.with_objects = [v if isinstance(v, ObjectId) else ObjectId(v) for v in self.with_objects]

        if not isinstance(self.provenances, list):
            self.provenances = [self.provenances] if self.provenances is not None else []
        self.provenances = [v if isinstance(v, ProvenanceInfo) else ProvenanceInfo(**as_dict(v)) for v in self.provenances]

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class Association(YAMLRoot):
    """
    An abstract grouping for different kinds of evidence-associated provenance
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["Association"]
    class_class_curie: ClassVar[str] = "gocam:Association"
    class_name: ClassVar[str] = "Association"
    class_model_uri: ClassVar[URIRef] = GOCAM.Association

    type: Optional[str] = None
    evidence: Optional[Union[Union[dict, EvidenceItem], list[Union[dict, EvidenceItem]]]] = empty_list()
    provenances: Optional[Union[Union[dict, "ProvenanceInfo"], list[Union[dict, "ProvenanceInfo"]]]] = empty_list()

    def __post_init__(self, *_: str, **kwargs: Any):
        self.type = str(self.class_name)

        if not isinstance(self.evidence, list):
            self.evidence = [self.evidence] if self.evidence is not None else []
        self.evidence = [v if isinstance(v, EvidenceItem) else EvidenceItem(**as_dict(v)) for v in self.evidence]

        if not isinstance(self.provenances, list):
            self.provenances = [self.provenances] if self.provenances is not None else []
        self.provenances = [v if isinstance(v, ProvenanceInfo) else ProvenanceInfo(**as_dict(v)) for v in self.provenances]

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


    def __new__(cls, *args, **kwargs):

        type_designator = "type"
        if not type_designator in kwargs:
            return super().__new__(cls,*args,**kwargs)
        else:
            type_designator_value = kwargs[type_designator]
            target_cls = cls._class_for("class_name", type_designator_value)


            if target_cls is None:
                raise ValueError(f"Wrong type designator value: class {cls.__name__} "
                                 f"has no subclass with ['class_name']='{kwargs[type_designator]}'")
            return super().__new__(target_cls,*args,**kwargs)



@dataclass(repr=False)
class EnabledByAssociation(Association):
    """
    An association between an activity and the gene product or complex or set of potential gene products
    that carry out that activity.
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["EnabledByAssociation"]
    class_class_curie: ClassVar[str] = "gocam:EnabledByAssociation"
    class_name: ClassVar[str] = "EnabledByAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.EnabledByAssociation

    term: Optional[Union[str, InformationBiomacromoleculeTermObjectId]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.term is not None and not isinstance(self.term, InformationBiomacromoleculeTermObjectId):
            self.term = InformationBiomacromoleculeTermObjectId(self.term)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass(repr=False)
class EnabledByGeneProductAssociation(EnabledByAssociation):
    """
    An association between an activity and an individual gene product
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["EnabledByGeneProductAssociation"]
    class_class_curie: ClassVar[str] = "gocam:EnabledByGeneProductAssociation"
    class_name: ClassVar[str] = "EnabledByGeneProductAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.EnabledByGeneProductAssociation

    term: Optional[Union[str, GeneProductTermObjectId]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.term is not None and not isinstance(self.term, GeneProductTermObjectId):
            self.term = GeneProductTermObjectId(self.term)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass(repr=False)
class EnabledByProteinComplexAssociation(EnabledByAssociation):
    """
    An association between an activity and a protein complex, where the complex carries out the activity. This should
    only be used when the activity cannot be attributed to an individual member of the complex, but instead the
    function is an emergent property of the complex.
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["EnabledByProteinComplexAssociation"]
    class_class_curie: ClassVar[str] = "gocam:EnabledByProteinComplexAssociation"
    class_name: ClassVar[str] = "EnabledByProteinComplexAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.EnabledByProteinComplexAssociation

    members: Optional[Union[Union[str, GeneProductTermObjectId], list[Union[str, GeneProductTermObjectId]]]] = empty_list()
    term: Optional[Union[str, ProteinComplexTermObjectId]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if not isinstance(self.members, list):
            self.members = [self.members] if self.members is not None else []
        self.members = [v if isinstance(v, GeneProductTermObjectId) else GeneProductTermObjectId(v) for v in self.members]

        if self.term is not None and not isinstance(self.term, ProteinComplexTermObjectId):
            self.term = ProteinComplexTermObjectId(self.term)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass(repr=False)
class CausalAssociation(Association):
    """
    A causal association between two activities
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["CausalAssociation"]
    class_class_curie: ClassVar[str] = "gocam:CausalAssociation"
    class_name: ClassVar[str] = "CausalAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.CausalAssociation

    predicate: Optional[Union[str, PredicateTermObjectId]] = None
    downstream_activity: Optional[Union[str, ActivityId]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.predicate is not None and not isinstance(self.predicate, PredicateTermObjectId):
            self.predicate = PredicateTermObjectId(self.predicate)

        if self.downstream_activity is not None and not isinstance(self.downstream_activity, ActivityId):
            self.downstream_activity = ActivityId(self.downstream_activity)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass(repr=False)
class TermAssociation(Association):
    """
    An association between an activity and a term, potentially with extensions. This is an abstract class for grouping
    purposes, it should not be directly instantiated, instead a subclass should be instantiated.
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["TermAssociation"]
    class_class_curie: ClassVar[str] = "gocam:TermAssociation"
    class_name: ClassVar[str] = "TermAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.TermAssociation

    term: Optional[Union[str, TermObjectId]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.term is not None and not isinstance(self.term, TermObjectId):
            self.term = TermObjectId(self.term)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass(repr=False)
class MolecularFunctionAssociation(TermAssociation):
    """
    An association between an activity and a molecular function term
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["MolecularFunctionAssociation"]
    class_class_curie: ClassVar[str] = "gocam:MolecularFunctionAssociation"
    class_name: ClassVar[str] = "MolecularFunctionAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.MolecularFunctionAssociation

    term: Optional[Union[str, MolecularFunctionTermObjectId]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.term is not None and not isinstance(self.term, MolecularFunctionTermObjectId):
            self.term = MolecularFunctionTermObjectId(self.term)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass(repr=False)
class BiologicalProcessAssociation(TermAssociation):
    """
    An association between an activity and a biological process term
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["BiologicalProcessAssociation"]
    class_class_curie: ClassVar[str] = "gocam:BiologicalProcessAssociation"
    class_name: ClassVar[str] = "BiologicalProcessAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.BiologicalProcessAssociation

    happens_during: Optional[Union[str, PhaseTermObjectId]] = None
    part_of: Optional[Union[dict, "BiologicalProcessAssociation"]] = None
    term: Optional[Union[str, BiologicalProcessTermObjectId]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.happens_during is not None and not isinstance(self.happens_during, PhaseTermObjectId):
            self.happens_during = PhaseTermObjectId(self.happens_during)

        if self.part_of is not None and not isinstance(self.part_of, BiologicalProcessAssociation):
            self.part_of = BiologicalProcessAssociation(**as_dict(self.part_of))

        if self.term is not None and not isinstance(self.term, BiologicalProcessTermObjectId):
            self.term = BiologicalProcessTermObjectId(self.term)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass(repr=False)
class CellularAnatomicalEntityAssociation(TermAssociation):
    """
    An association between an activity and a cellular anatomical entity term
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["CellularAnatomicalEntityAssociation"]
    class_class_curie: ClassVar[str] = "gocam:CellularAnatomicalEntityAssociation"
    class_name: ClassVar[str] = "CellularAnatomicalEntityAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.CellularAnatomicalEntityAssociation

    part_of: Optional[Union[dict, "CellTypeAssociation"]] = None
    term: Optional[Union[str, CellularAnatomicalEntityTermObjectId]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.part_of is not None and not isinstance(self.part_of, CellTypeAssociation):
            self.part_of = CellTypeAssociation(**as_dict(self.part_of))

        if self.term is not None and not isinstance(self.term, CellularAnatomicalEntityTermObjectId):
            self.term = CellularAnatomicalEntityTermObjectId(self.term)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass(repr=False)
class CellTypeAssociation(TermAssociation):
    """
    An association between an activity and a cell type term
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["CellTypeAssociation"]
    class_class_curie: ClassVar[str] = "gocam:CellTypeAssociation"
    class_name: ClassVar[str] = "CellTypeAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.CellTypeAssociation

    part_of: Optional[Union[dict, "GrossAnatomyAssociation"]] = None
    term: Optional[Union[str, CellTypeTermObjectId]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.part_of is not None and not isinstance(self.part_of, GrossAnatomyAssociation):
            self.part_of = GrossAnatomyAssociation(**as_dict(self.part_of))

        if self.term is not None and not isinstance(self.term, CellTypeTermObjectId):
            self.term = CellTypeTermObjectId(self.term)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass(repr=False)
class GrossAnatomyAssociation(TermAssociation):
    """
    An association between an activity and a gross anatomical structure term
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["GrossAnatomyAssociation"]
    class_class_curie: ClassVar[str] = "gocam:GrossAnatomyAssociation"
    class_name: ClassVar[str] = "GrossAnatomyAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.GrossAnatomyAssociation

    part_of: Optional[Union[dict, "GrossAnatomyAssociation"]] = None
    term: Optional[Union[str, GrossAnatomicalStructureTermObjectId]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.part_of is not None and not isinstance(self.part_of, GrossAnatomyAssociation):
            self.part_of = GrossAnatomyAssociation(**as_dict(self.part_of))

        if self.term is not None and not isinstance(self.term, GrossAnatomicalStructureTermObjectId):
            self.term = GrossAnatomicalStructureTermObjectId(self.term)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass(repr=False)
class MoleculeAssociation(TermAssociation):
    """
    An association between an activity and a molecule term
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["MoleculeAssociation"]
    class_class_curie: ClassVar[str] = "gocam:MoleculeAssociation"
    class_name: ClassVar[str] = "MoleculeAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.MoleculeAssociation

    term: Optional[Union[str, MoleculeTermObjectId]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.term is not None and not isinstance(self.term, MoleculeTermObjectId):
            self.term = MoleculeTermObjectId(self.term)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass(repr=False)
class Object(YAMLRoot):
    """
    An abstract class for all identified objects in a model
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["Object"]
    class_class_curie: ClassVar[str] = "gocam:Object"
    class_name: ClassVar[str] = "Object"
    class_model_uri: ClassVar[URIRef] = GOCAM.Object

    id: Union[str, ObjectId] = None
    label: Optional[str] = None
    type: Optional[Union[str, URIorCURIE]] = None
    obsolete: Optional[Union[bool, Bool]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, ObjectId):
            self.id = ObjectId(self.id)

        if self.label is not None and not isinstance(self.label, str):
            self.label = str(self.label)

        self.type = str(self.class_class_curie)

        if self.obsolete is not None and not isinstance(self.obsolete, Bool):
            self.obsolete = Bool(self.obsolete)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


    def __new__(cls, *args, **kwargs):

        type_designator = "type"
        if not type_designator in kwargs:
            return super().__new__(cls,*args,**kwargs)
        else:
            type_designator_value = kwargs[type_designator]
            target_cls = cls._class_for("class_class_curie", type_designator_value)


            if target_cls is None:
                target_cls = cls._class_for("class_class_uri", type_designator_value)


            if target_cls is None:
                target_cls = cls._class_for("class_model_uri", type_designator_value)


            if target_cls is None:
                raise ValueError(f"Wrong type designator value: class {cls.__name__} "
                                 f"has no subclass with ['class_class_curie', 'class_class_uri', 'class_model_uri']='{kwargs[type_designator]}'")
            return super().__new__(target_cls,*args,**kwargs)



@dataclass(repr=False)
class TermObject(Object):
    """
    An abstract class for all ontology term objects
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["TermObject"]
    class_class_curie: ClassVar[str] = "gocam:TermObject"
    class_name: ClassVar[str] = "TermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.TermObject

    id: Union[str, TermObjectId] = None

    def __post_init__(self, *_: str, **kwargs: Any):

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass(repr=False)
class PublicationObject(Object):
    """
    An object that represents a publication or other kind of reference
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["PublicationObject"]
    class_class_curie: ClassVar[str] = "gocam:PublicationObject"
    class_name: ClassVar[str] = "PublicationObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.PublicationObject

    id: Union[str, PublicationObjectId] = None
    abstract_text: Optional[str] = None
    full_text: Optional[str] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, PublicationObjectId):
            self.id = PublicationObjectId(self.id)

        if self.abstract_text is not None and not isinstance(self.abstract_text, str):
            self.abstract_text = str(self.abstract_text)

        if self.full_text is not None and not isinstance(self.full_text, str):
            self.full_text = str(self.full_text)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass(repr=False)
class EvidenceTermObject(TermObject):
    """
    A term object that represents an evidence term from ECO. Only ECO terms that map up to a GO GAF evidence code
    should be used.
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["EvidenceTermObject"]
    class_class_curie: ClassVar[str] = "gocam:EvidenceTermObject"
    class_name: ClassVar[str] = "EvidenceTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.EvidenceTermObject

    id: Union[str, EvidenceTermObjectId] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, EvidenceTermObjectId):
            self.id = EvidenceTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass(repr=False)
class MolecularFunctionTermObject(TermObject):
    """
    A term object that represents a molecular function term from GO
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["MolecularFunctionTermObject"]
    class_class_curie: ClassVar[str] = "gocam:MolecularFunctionTermObject"
    class_name: ClassVar[str] = "MolecularFunctionTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.MolecularFunctionTermObject

    id: Union[str, MolecularFunctionTermObjectId] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, MolecularFunctionTermObjectId):
            self.id = MolecularFunctionTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass(repr=False)
class BiologicalProcessTermObject(TermObject):
    """
    A term object that represents a biological process term from GO
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["BiologicalProcessTermObject"]
    class_class_curie: ClassVar[str] = "gocam:BiologicalProcessTermObject"
    class_name: ClassVar[str] = "BiologicalProcessTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.BiologicalProcessTermObject

    id: Union[str, BiologicalProcessTermObjectId] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, BiologicalProcessTermObjectId):
            self.id = BiologicalProcessTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass(repr=False)
class CellularAnatomicalEntityTermObject(TermObject):
    """
    A term object that represents a cellular anatomical entity term from GO
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["CellularAnatomicalEntityTermObject"]
    class_class_curie: ClassVar[str] = "gocam:CellularAnatomicalEntityTermObject"
    class_name: ClassVar[str] = "CellularAnatomicalEntityTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.CellularAnatomicalEntityTermObject

    id: Union[str, CellularAnatomicalEntityTermObjectId] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, CellularAnatomicalEntityTermObjectId):
            self.id = CellularAnatomicalEntityTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass(repr=False)
class MoleculeTermObject(TermObject):
    """
    A term object that represents a molecule term from CHEBI or UniProtKB
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["MoleculeTermObject"]
    class_class_curie: ClassVar[str] = "gocam:MoleculeTermObject"
    class_name: ClassVar[str] = "MoleculeTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.MoleculeTermObject

    id: Union[str, MoleculeTermObjectId] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, MoleculeTermObjectId):
            self.id = MoleculeTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass(repr=False)
class CellTypeTermObject(TermObject):
    """
    A term object that represents a cell type term from CL
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["CellTypeTermObject"]
    class_class_curie: ClassVar[str] = "gocam:CellTypeTermObject"
    class_name: ClassVar[str] = "CellTypeTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.CellTypeTermObject

    id: Union[str, CellTypeTermObjectId] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, CellTypeTermObjectId):
            self.id = CellTypeTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass(repr=False)
class GrossAnatomicalStructureTermObject(TermObject):
    """
    A term object that represents a gross anatomical structure term from UBERON
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["GrossAnatomicalStructureTermObject"]
    class_class_curie: ClassVar[str] = "gocam:GrossAnatomicalStructureTermObject"
    class_name: ClassVar[str] = "GrossAnatomicalStructureTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.GrossAnatomicalStructureTermObject

    id: Union[str, GrossAnatomicalStructureTermObjectId] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, GrossAnatomicalStructureTermObjectId):
            self.id = GrossAnatomicalStructureTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass(repr=False)
class PhaseTermObject(TermObject):
    """
    A term object that represents a phase term from GO or UBERON
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["PhaseTermObject"]
    class_class_curie: ClassVar[str] = "gocam:PhaseTermObject"
    class_name: ClassVar[str] = "PhaseTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.PhaseTermObject

    id: Union[str, PhaseTermObjectId] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, PhaseTermObjectId):
            self.id = PhaseTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass(repr=False)
class InformationBiomacromoleculeTermObject(TermObject):
    """
    An abstract class for all information biomacromolecule term objects
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["InformationBiomacromoleculeTermObject"]
    class_class_curie: ClassVar[str] = "gocam:InformationBiomacromoleculeTermObject"
    class_name: ClassVar[str] = "InformationBiomacromoleculeTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.InformationBiomacromoleculeTermObject

    id: Union[str, InformationBiomacromoleculeTermObjectId] = None

    def __post_init__(self, *_: str, **kwargs: Any):

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass(repr=False)
class GeneProductTermObject(InformationBiomacromoleculeTermObject):
    """
    A term object that represents a gene product term from GO or UniProtKB
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["GeneProductTermObject"]
    class_class_curie: ClassVar[str] = "gocam:GeneProductTermObject"
    class_name: ClassVar[str] = "GeneProductTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.GeneProductTermObject

    id: Union[str, GeneProductTermObjectId] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, GeneProductTermObjectId):
            self.id = GeneProductTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass(repr=False)
class ProteinComplexTermObject(InformationBiomacromoleculeTermObject):
    """
    A term object that represents a protein complex term from GO
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["ProteinComplexTermObject"]
    class_class_curie: ClassVar[str] = "gocam:ProteinComplexTermObject"
    class_name: ClassVar[str] = "ProteinComplexTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.ProteinComplexTermObject

    id: Union[str, ProteinComplexTermObjectId] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, ProteinComplexTermObjectId):
            self.id = ProteinComplexTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass(repr=False)
class TaxonTermObject(TermObject):
    """
    A term object that represents a taxon term from NCBITaxon
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["TaxonTermObject"]
    class_class_curie: ClassVar[str] = "gocam:TaxonTermObject"
    class_name: ClassVar[str] = "TaxonTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.TaxonTermObject

    id: Union[str, TaxonTermObjectId] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, TaxonTermObjectId):
            self.id = TaxonTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass(repr=False)
class PredicateTermObject(TermObject):
    """
    A term object that represents a taxon term from NCBITaxon
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["PredicateTermObject"]
    class_class_curie: ClassVar[str] = "gocam:PredicateTermObject"
    class_name: ClassVar[str] = "PredicateTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.PredicateTermObject

    id: Union[str, PredicateTermObjectId] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, PredicateTermObjectId):
            self.id = PredicateTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass(repr=False)
class ProvenanceInfo(YAMLRoot):
    """
    Provenance information for an object
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["ProvenanceInfo"]
    class_class_curie: ClassVar[str] = "gocam:ProvenanceInfo"
    class_name: ClassVar[str] = "ProvenanceInfo"
    class_model_uri: ClassVar[URIRef] = GOCAM.ProvenanceInfo

    contributor: Optional[Union[str, list[str]]] = empty_list()
    created: Optional[str] = None
    date: Optional[str] = None
    provided_by: Optional[Union[str, list[str]]] = empty_list()

    def __post_init__(self, *_: str, **kwargs: Any):
        if not isinstance(self.contributor, list):
            self.contributor = [self.contributor] if self.contributor is not None else []
        self.contributor = [v if isinstance(v, str) else str(v) for v in self.contributor]

        if self.created is not None and not isinstance(self.created, str):
            self.created = str(self.created)

        if self.date is not None and not isinstance(self.date, str):
            self.date = str(self.date)

        if not isinstance(self.provided_by, list):
            self.provided_by = [self.provided_by] if self.provided_by is not None else []
        self.provided_by = [v if isinstance(v, str) else str(v) for v in self.provided_by]

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class QueryIndex(YAMLRoot):
    """
    An index that is optionally placed on a model in order to support common query or index operations. Note that this
    index is not typically populated in the working transactional store for a model, it is derived via computation
    from core primary model information.
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["QueryIndex"]
    class_class_curie: ClassVar[str] = "gocam:QueryIndex"
    class_name: ClassVar[str] = "QueryIndex"
    class_model_uri: ClassVar[URIRef] = GOCAM.QueryIndex

    number_of_activities: Optional[int] = None
    number_of_enabled_by_terms: Optional[int] = None
    number_of_causal_associations: Optional[int] = None
    length_of_longest_causal_association_path: Optional[int] = None
    number_of_strongly_connected_components: Optional[int] = None
    flattened_references: Optional[Union[dict[Union[str, PublicationObjectId], Union[dict, PublicationObject]], list[Union[dict, PublicationObject]]]] = empty_dict()
    model_activity_molecular_function_terms: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_activity_molecular_function_closure: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_activity_molecular_function_rollup: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_activity_occurs_in_terms: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_activity_occurs_in_closure: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_activity_occurs_in_rollup: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_activity_enabled_by_terms: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_activity_enabled_by_closure: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_activity_enabled_by_rollup: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_activity_part_of_terms: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_activity_part_of_closure: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_activity_part_of_rollup: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_activity_has_input_terms: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_activity_has_input_closure: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_activity_has_input_rollup: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_taxon: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_taxon_closure: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    model_taxon_rollup: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    annoton_terms: Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]] = empty_dict()
    start_activities: Optional[Union[Union[str, ActivityId], list[Union[str, ActivityId]]]] = empty_list()
    end_activities: Optional[Union[Union[str, ActivityId], list[Union[str, ActivityId]]]] = empty_list()
    intermediate_activities: Optional[Union[Union[str, ActivityId], list[Union[str, ActivityId]]]] = empty_list()
    singleton_activities: Optional[Union[Union[str, ActivityId], list[Union[str, ActivityId]]]] = empty_list()
    number_of_start_activities: Optional[int] = None
    number_of_end_activities: Optional[int] = None
    number_of_intermediate_activities: Optional[int] = None
    number_of_singleton_activities: Optional[int] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.number_of_activities is not None and not isinstance(self.number_of_activities, int):
            self.number_of_activities = int(self.number_of_activities)

        if self.number_of_enabled_by_terms is not None and not isinstance(self.number_of_enabled_by_terms, int):
            self.number_of_enabled_by_terms = int(self.number_of_enabled_by_terms)

        if self.number_of_causal_associations is not None and not isinstance(self.number_of_causal_associations, int):
            self.number_of_causal_associations = int(self.number_of_causal_associations)

        if self.length_of_longest_causal_association_path is not None and not isinstance(self.length_of_longest_causal_association_path, int):
            self.length_of_longest_causal_association_path = int(self.length_of_longest_causal_association_path)

        if self.number_of_strongly_connected_components is not None and not isinstance(self.number_of_strongly_connected_components, int):
            self.number_of_strongly_connected_components = int(self.number_of_strongly_connected_components)

        self._normalize_inlined_as_list(slot_name="flattened_references", slot_type=PublicationObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_activity_molecular_function_terms", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_activity_molecular_function_closure", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_activity_molecular_function_rollup", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_activity_occurs_in_terms", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_activity_occurs_in_closure", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_activity_occurs_in_rollup", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_activity_enabled_by_terms", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_activity_enabled_by_closure", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_activity_enabled_by_rollup", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_activity_part_of_terms", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_activity_part_of_closure", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_activity_part_of_rollup", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_activity_has_input_terms", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_activity_has_input_closure", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_activity_has_input_rollup", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_taxon", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_taxon_closure", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="model_taxon_rollup", slot_type=TermObject, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="annoton_terms", slot_type=TermObject, key_name="id", keyed=True)

        if not isinstance(self.start_activities, list):
            self.start_activities = [self.start_activities] if self.start_activities is not None else []
        self.start_activities = [v if isinstance(v, ActivityId) else ActivityId(v) for v in self.start_activities]

        if not isinstance(self.end_activities, list):
            self.end_activities = [self.end_activities] if self.end_activities is not None else []
        self.end_activities = [v if isinstance(v, ActivityId) else ActivityId(v) for v in self.end_activities]

        if not isinstance(self.intermediate_activities, list):
            self.intermediate_activities = [self.intermediate_activities] if self.intermediate_activities is not None else []
        self.intermediate_activities = [v if isinstance(v, ActivityId) else ActivityId(v) for v in self.intermediate_activities]

        if not isinstance(self.singleton_activities, list):
            self.singleton_activities = [self.singleton_activities] if self.singleton_activities is not None else []
        self.singleton_activities = [v if isinstance(v, ActivityId) else ActivityId(v) for v in self.singleton_activities]

        if self.number_of_start_activities is not None and not isinstance(self.number_of_start_activities, int):
            self.number_of_start_activities = int(self.number_of_start_activities)

        if self.number_of_end_activities is not None and not isinstance(self.number_of_end_activities, int):
            self.number_of_end_activities = int(self.number_of_end_activities)

        if self.number_of_intermediate_activities is not None and not isinstance(self.number_of_intermediate_activities, int):
            self.number_of_intermediate_activities = int(self.number_of_intermediate_activities)

        if self.number_of_singleton_activities is not None and not isinstance(self.number_of_singleton_activities, int):
            self.number_of_singleton_activities = int(self.number_of_singleton_activities)

        super().__post_init__(**kwargs)


# Enumerations
class ModelStateEnum(EnumDefinitionImpl):
    """
    A term describing where the model is in the development life cycle.
    """
    development = PermissibleValue(
        text="development",
        description="""Used when the curator is still working on the model. Edits are still being made, and the information in the model is not yet guaranteed to be accurate or complete. The model should not be displayed in end-user facing websites, unless it is made clear that the model is a work in progress.""")
    production = PermissibleValue(
        text="production",
        description="""Used when the curator has declared the model is ready for public consumption. Edits might still be performed on the model in future, but the information in the model is believed to be both accurate and reasonably complete. The model may be displayed in public websites.""")
    delete = PermissibleValue(
        text="delete",
        description="When the curator has marked for future deletion.")
    review = PermissibleValue(
        text="review",
        description="The model has been marked for curator review.")
    internal_test = PermissibleValue(
        text="internal_test",
        description="The model is not intended for use public use; it is likely to be used for internal testing.")
    closed = PermissibleValue(
        text="closed",
        description="TBD")

    _defn = EnumDefinition(
        name="ModelStateEnum",
        description="A term describing where the model is in the development life cycle.",
    )

class InformationBiomacromoleculeCategory(EnumDefinitionImpl):
    """
    A term describing the type of the enabler of an activity.
    """
    GeneOrReferenceProtein = PermissibleValue(
        text="GeneOrReferenceProtein",
        meaning=GOCAM["biolink.GeneOrGeneProduct"])
    ProteinIsoform = PermissibleValue(text="ProteinIsoform")
    MacromolecularComplex = PermissibleValue(text="MacromolecularComplex")
    Unknown = PermissibleValue(text="Unknown")

    _defn = EnumDefinition(
        name="InformationBiomacromoleculeCategory",
        description="A term describing the type of the enabler of an activity.",
    )

class CausalPredicateEnum(EnumDefinitionImpl):
    """
    A term describing the causal relationship between two activities. All terms are drawn from the "causally upstream
    or within" (RO:0002418) branch of the Relation Ontology (RO).
    """
    regulates = PermissibleValue(
        text="regulates",
        meaning=RO["0002211"])

    _defn = EnumDefinition(
        name="CausalPredicateEnum",
        description="""A term describing the causal relationship between two activities. All terms are drawn from the \"causally upstream or within\" (RO:0002418) branch of the Relation Ontology (RO).""",
    )

    @classmethod
    def _addvals(cls):
        setattr(cls, "causally upstream of, positive effect",
            PermissibleValue(
                text="causally upstream of, positive effect",
                meaning=RO["0002304"]))
        setattr(cls, "causally upstream of, negative effect",
            PermissibleValue(
                text="causally upstream of, negative effect",
                meaning=RO["0002305"]))
        setattr(cls, "causally upstream of",
            PermissibleValue(
                text="causally upstream of",
                meaning=RO["0002411"]))
        setattr(cls, "immediately causally upstream of",
            PermissibleValue(
                text="immediately causally upstream of",
                meaning=RO["0002412"]))
        setattr(cls, "causally upstream of or within",
            PermissibleValue(
                text="causally upstream of or within",
                meaning=RO["0002418"]))
        setattr(cls, "causally upstream of or within, negative effect",
            PermissibleValue(
                text="causally upstream of or within, negative effect",
                meaning=RO["0004046"]))
        setattr(cls, "causally upstream of or within, positive effect",
            PermissibleValue(
                text="causally upstream of or within, positive effect",
                meaning=RO["0004047"]))
        setattr(cls, "negatively regulates",
            PermissibleValue(
                text="negatively regulates",
                meaning=RO["0002212"]))
        setattr(cls, "positively regulates",
            PermissibleValue(
                text="positively regulates",
                meaning=RO["0002213"]))
        setattr(cls, "provides input for",
            PermissibleValue(
                text="provides input for",
                meaning=RO["0002413"]))
        setattr(cls, "removes input for",
            PermissibleValue(
                text="removes input for",
                meaning=RO["0012010"]))

class EvidenceCodeEnum(EnumDefinitionImpl):
    """
    A term from the subset of ECO that maps up to a GAF evidence code
    """
    _defn = EnumDefinition(
        name="EvidenceCodeEnum",
        description="A term from the subset of ECO that maps up to a GAF evidence code",
    )

class CellularAnatomicalEntityEnum(EnumDefinitionImpl):
    """
    A term from the subset of the cellular anatomical entity branch of GO CC
    """
    _defn = EnumDefinition(
        name="CellularAnatomicalEntityEnum",
        description="A term from the subset of the cellular anatomical entity branch of GO CC",
    )

class PhaseEnum(EnumDefinitionImpl):
    """
    A term from either the phase branch of GO or the phase branch of an anatomy ontology
    """
    _defn = EnumDefinition(
        name="PhaseEnum",
        description="A term from either the phase branch of GO or the phase branch of an anatomy ontology",
    )

# Slots
class slots:
    pass

slots.model__id = Slot(uri=GOCAM.id, name="model__id", curie=GOCAM.curie('id'),
                   model_uri=GOCAM.model__id, domain=None, range=URIRef)

slots.model__title = Slot(uri=DCT.title, name="model__title", curie=DCT.curie('title'),
                   model_uri=GOCAM.model__title, domain=None, range=str)

slots.model__taxon = Slot(uri=GOCAM.taxon, name="model__taxon", curie=GOCAM.curie('taxon'),
                   model_uri=GOCAM.model__taxon, domain=None, range=Optional[Union[str, TaxonTermObjectId]])

slots.model__additional_taxa = Slot(uri=GOCAM.additional_taxa, name="model__additional_taxa", curie=GOCAM.curie('additional_taxa'),
                   model_uri=GOCAM.model__additional_taxa, domain=None, range=Optional[Union[Union[str, TaxonTermObjectId], list[Union[str, TaxonTermObjectId]]]])

slots.model__status = Slot(uri=PAV.status, name="model__status", curie=PAV.curie('status'),
                   model_uri=GOCAM.model__status, domain=None, range=Optional[Union[str, "ModelStateEnum"]])

slots.model__comments = Slot(uri=RDFS.comment, name="model__comments", curie=RDFS.curie('comment'),
                   model_uri=GOCAM.model__comments, domain=None, range=Optional[Union[str, list[str]]])

slots.model__activities = Slot(uri=GOCAM.activities, name="model__activities", curie=GOCAM.curie('activities'),
                   model_uri=GOCAM.model__activities, domain=None, range=Optional[Union[dict[Union[str, ActivityId], Union[dict, Activity]], list[Union[dict, Activity]]]])

slots.model__objects = Slot(uri=GOCAM.objects, name="model__objects", curie=GOCAM.curie('objects'),
                   model_uri=GOCAM.model__objects, domain=None, range=Optional[Union[dict[Union[str, ObjectId], Union[dict, Object]], list[Union[dict, Object]]]])

slots.model__provenances = Slot(uri=GOCAM.provenances, name="model__provenances", curie=GOCAM.curie('provenances'),
                   model_uri=GOCAM.model__provenances, domain=None, range=Optional[Union[Union[dict, ProvenanceInfo], list[Union[dict, ProvenanceInfo]]]])

slots.model__query_index = Slot(uri=GOCAM.query_index, name="model__query_index", curie=GOCAM.curie('query_index'),
                   model_uri=GOCAM.model__query_index, domain=None, range=Optional[Union[dict, QueryIndex]])

slots.activity__id = Slot(uri=GOCAM.id, name="activity__id", curie=GOCAM.curie('id'),
                   model_uri=GOCAM.activity__id, domain=None, range=URIRef)

slots.activity__enabled_by = Slot(uri=GOCAM.enabled_by, name="activity__enabled_by", curie=GOCAM.curie('enabled_by'),
                   model_uri=GOCAM.activity__enabled_by, domain=None, range=Optional[Union[dict, EnabledByAssociation]])

slots.activity__molecular_function = Slot(uri=GOCAM.molecular_function, name="activity__molecular_function", curie=GOCAM.curie('molecular_function'),
                   model_uri=GOCAM.activity__molecular_function, domain=None, range=Optional[Union[dict, MolecularFunctionAssociation]])

slots.activity__occurs_in = Slot(uri=GOCAM.occurs_in, name="activity__occurs_in", curie=GOCAM.curie('occurs_in'),
                   model_uri=GOCAM.activity__occurs_in, domain=None, range=Optional[Union[dict, CellularAnatomicalEntityAssociation]])

slots.activity__part_of = Slot(uri=GOCAM.part_of, name="activity__part_of", curie=GOCAM.curie('part_of'),
                   model_uri=GOCAM.activity__part_of, domain=None, range=Optional[Union[dict, BiologicalProcessAssociation]])

slots.activity__has_input = Slot(uri=GOCAM.has_input, name="activity__has_input", curie=GOCAM.curie('has_input'),
                   model_uri=GOCAM.activity__has_input, domain=None, range=Optional[Union[Union[dict, MoleculeAssociation], list[Union[dict, MoleculeAssociation]]]])

slots.activity__has_primary_input = Slot(uri=GOCAM.has_primary_input, name="activity__has_primary_input", curie=GOCAM.curie('has_primary_input'),
                   model_uri=GOCAM.activity__has_primary_input, domain=None, range=Optional[Union[dict, MoleculeAssociation]])

slots.activity__has_output = Slot(uri=GOCAM.has_output, name="activity__has_output", curie=GOCAM.curie('has_output'),
                   model_uri=GOCAM.activity__has_output, domain=None, range=Optional[Union[Union[dict, MoleculeAssociation], list[Union[dict, MoleculeAssociation]]]])

slots.activity__has_primary_output = Slot(uri=GOCAM.has_primary_output, name="activity__has_primary_output", curie=GOCAM.curie('has_primary_output'),
                   model_uri=GOCAM.activity__has_primary_output, domain=None, range=Optional[Union[dict, MoleculeAssociation]])

slots.activity__causal_associations = Slot(uri=GOCAM.causal_associations, name="activity__causal_associations", curie=GOCAM.curie('causal_associations'),
                   model_uri=GOCAM.activity__causal_associations, domain=None, range=Optional[Union[Union[dict, CausalAssociation], list[Union[dict, CausalAssociation]]]])

slots.activity__provenances = Slot(uri=GOCAM.provenances, name="activity__provenances", curie=GOCAM.curie('provenances'),
                   model_uri=GOCAM.activity__provenances, domain=None, range=Optional[Union[Union[dict, ProvenanceInfo], list[Union[dict, ProvenanceInfo]]]])

slots.evidenceItem__term = Slot(uri=GOCAM.term, name="evidenceItem__term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.evidenceItem__term, domain=None, range=Optional[Union[str, EvidenceTermObjectId]])

slots.evidenceItem__reference = Slot(uri=GOCAM.reference, name="evidenceItem__reference", curie=GOCAM.curie('reference'),
                   model_uri=GOCAM.evidenceItem__reference, domain=None, range=Optional[Union[str, PublicationObjectId]])

slots.evidenceItem__with_objects = Slot(uri=GOCAM.with_objects, name="evidenceItem__with_objects", curie=GOCAM.curie('with_objects'),
                   model_uri=GOCAM.evidenceItem__with_objects, domain=None, range=Optional[Union[Union[str, ObjectId], list[Union[str, ObjectId]]]])

slots.evidenceItem__provenances = Slot(uri=GOCAM.provenances, name="evidenceItem__provenances", curie=GOCAM.curie('provenances'),
                   model_uri=GOCAM.evidenceItem__provenances, domain=None, range=Optional[Union[Union[dict, ProvenanceInfo], list[Union[dict, ProvenanceInfo]]]])

slots.association__type = Slot(uri=GOCAM.type, name="association__type", curie=GOCAM.curie('type'),
                   model_uri=GOCAM.association__type, domain=None, range=Optional[str])

slots.association__evidence = Slot(uri=GOCAM.evidence, name="association__evidence", curie=GOCAM.curie('evidence'),
                   model_uri=GOCAM.association__evidence, domain=None, range=Optional[Union[Union[dict, EvidenceItem], list[Union[dict, EvidenceItem]]]])

slots.association__provenances = Slot(uri=GOCAM.provenances, name="association__provenances", curie=GOCAM.curie('provenances'),
                   model_uri=GOCAM.association__provenances, domain=None, range=Optional[Union[Union[dict, ProvenanceInfo], list[Union[dict, ProvenanceInfo]]]])

slots.enabledByAssociation__term = Slot(uri=GOCAM.term, name="enabledByAssociation__term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.enabledByAssociation__term, domain=None, range=Optional[Union[str, InformationBiomacromoleculeTermObjectId]])

slots.enabledByProteinComplexAssociation__members = Slot(uri=GOCAM.members, name="enabledByProteinComplexAssociation__members", curie=GOCAM.curie('members'),
                   model_uri=GOCAM.enabledByProteinComplexAssociation__members, domain=None, range=Optional[Union[Union[str, GeneProductTermObjectId], list[Union[str, GeneProductTermObjectId]]]])

slots.causalAssociation__predicate = Slot(uri=GOCAM.predicate, name="causalAssociation__predicate", curie=GOCAM.curie('predicate'),
                   model_uri=GOCAM.causalAssociation__predicate, domain=None, range=Optional[Union[str, PredicateTermObjectId]])

slots.causalAssociation__downstream_activity = Slot(uri=GOCAM.downstream_activity, name="causalAssociation__downstream_activity", curie=GOCAM.curie('downstream_activity'),
                   model_uri=GOCAM.causalAssociation__downstream_activity, domain=None, range=Optional[Union[str, ActivityId]])

slots.termAssociation__term = Slot(uri=GOCAM.term, name="termAssociation__term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.termAssociation__term, domain=None, range=Optional[Union[str, TermObjectId]])

slots.biologicalProcessAssociation__happens_during = Slot(uri=GOCAM.happens_during, name="biologicalProcessAssociation__happens_during", curie=GOCAM.curie('happens_during'),
                   model_uri=GOCAM.biologicalProcessAssociation__happens_during, domain=None, range=Optional[Union[str, PhaseTermObjectId]])

slots.biologicalProcessAssociation__part_of = Slot(uri=GOCAM.part_of, name="biologicalProcessAssociation__part_of", curie=GOCAM.curie('part_of'),
                   model_uri=GOCAM.biologicalProcessAssociation__part_of, domain=None, range=Optional[Union[dict, BiologicalProcessAssociation]])

slots.cellularAnatomicalEntityAssociation__part_of = Slot(uri=GOCAM.part_of, name="cellularAnatomicalEntityAssociation__part_of", curie=GOCAM.curie('part_of'),
                   model_uri=GOCAM.cellularAnatomicalEntityAssociation__part_of, domain=None, range=Optional[Union[dict, CellTypeAssociation]])

slots.cellTypeAssociation__part_of = Slot(uri=GOCAM.part_of, name="cellTypeAssociation__part_of", curie=GOCAM.curie('part_of'),
                   model_uri=GOCAM.cellTypeAssociation__part_of, domain=None, range=Optional[Union[dict, GrossAnatomyAssociation]])

slots.grossAnatomyAssociation__part_of = Slot(uri=GOCAM.part_of, name="grossAnatomyAssociation__part_of", curie=GOCAM.curie('part_of'),
                   model_uri=GOCAM.grossAnatomyAssociation__part_of, domain=None, range=Optional[Union[dict, GrossAnatomyAssociation]])

slots.object__id = Slot(uri=GOCAM.id, name="object__id", curie=GOCAM.curie('id'),
                   model_uri=GOCAM.object__id, domain=None, range=URIRef)

slots.object__label = Slot(uri=RDFS.label, name="object__label", curie=RDFS.curie('label'),
                   model_uri=GOCAM.object__label, domain=None, range=Optional[str])

slots.object__type = Slot(uri=GOCAM.type, name="object__type", curie=GOCAM.curie('type'),
                   model_uri=GOCAM.object__type, domain=None, range=Optional[Union[str, URIorCURIE]])

slots.object__obsolete = Slot(uri=GOCAM.obsolete, name="object__obsolete", curie=GOCAM.curie('obsolete'),
                   model_uri=GOCAM.object__obsolete, domain=None, range=Optional[Union[bool, Bool]])

slots.publicationObject__abstract_text = Slot(uri=GOCAM.abstract_text, name="publicationObject__abstract_text", curie=GOCAM.curie('abstract_text'),
                   model_uri=GOCAM.publicationObject__abstract_text, domain=None, range=Optional[str])

slots.publicationObject__full_text = Slot(uri=GOCAM.full_text, name="publicationObject__full_text", curie=GOCAM.curie('full_text'),
                   model_uri=GOCAM.publicationObject__full_text, domain=None, range=Optional[str])

slots.provenanceInfo__contributor = Slot(uri=DCT.contributor, name="provenanceInfo__contributor", curie=DCT.curie('contributor'),
                   model_uri=GOCAM.provenanceInfo__contributor, domain=None, range=Optional[Union[str, list[str]]])

slots.provenanceInfo__created = Slot(uri=DCT.created, name="provenanceInfo__created", curie=DCT.curie('created'),
                   model_uri=GOCAM.provenanceInfo__created, domain=None, range=Optional[str])

slots.provenanceInfo__date = Slot(uri=DCT.date, name="provenanceInfo__date", curie=DCT.curie('date'),
                   model_uri=GOCAM.provenanceInfo__date, domain=None, range=Optional[str])

slots.provenanceInfo__provided_by = Slot(uri=PAV.providedBy, name="provenanceInfo__provided_by", curie=PAV.curie('providedBy'),
                   model_uri=GOCAM.provenanceInfo__provided_by, domain=None, range=Optional[Union[str, list[str]]])

slots.queryIndex__number_of_activities = Slot(uri=GOCAM.number_of_activities, name="queryIndex__number_of_activities", curie=GOCAM.curie('number_of_activities'),
                   model_uri=GOCAM.queryIndex__number_of_activities, domain=None, range=Optional[int])

slots.queryIndex__number_of_enabled_by_terms = Slot(uri=GOCAM.number_of_enabled_by_terms, name="queryIndex__number_of_enabled_by_terms", curie=GOCAM.curie('number_of_enabled_by_terms'),
                   model_uri=GOCAM.queryIndex__number_of_enabled_by_terms, domain=None, range=Optional[int])

slots.queryIndex__number_of_causal_associations = Slot(uri=GOCAM.number_of_causal_associations, name="queryIndex__number_of_causal_associations", curie=GOCAM.curie('number_of_causal_associations'),
                   model_uri=GOCAM.queryIndex__number_of_causal_associations, domain=None, range=Optional[int])

slots.queryIndex__length_of_longest_causal_association_path = Slot(uri=GOCAM.length_of_longest_causal_association_path, name="queryIndex__length_of_longest_causal_association_path", curie=GOCAM.curie('length_of_longest_causal_association_path'),
                   model_uri=GOCAM.queryIndex__length_of_longest_causal_association_path, domain=None, range=Optional[int])

slots.queryIndex__number_of_strongly_connected_components = Slot(uri=GOCAM.number_of_strongly_connected_components, name="queryIndex__number_of_strongly_connected_components", curie=GOCAM.curie('number_of_strongly_connected_components'),
                   model_uri=GOCAM.queryIndex__number_of_strongly_connected_components, domain=None, range=Optional[int])

slots.queryIndex__flattened_references = Slot(uri=GOCAM.flattened_references, name="queryIndex__flattened_references", curie=GOCAM.curie('flattened_references'),
                   model_uri=GOCAM.queryIndex__flattened_references, domain=None, range=Optional[Union[dict[Union[str, PublicationObjectId], Union[dict, PublicationObject]], list[Union[dict, PublicationObject]]]])

slots.queryIndex__model_activity_molecular_function_terms = Slot(uri=GOCAM.model_activity_molecular_function_terms, name="queryIndex__model_activity_molecular_function_terms", curie=GOCAM.curie('model_activity_molecular_function_terms'),
                   model_uri=GOCAM.queryIndex__model_activity_molecular_function_terms, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_activity_molecular_function_closure = Slot(uri=GOCAM.model_activity_molecular_function_closure, name="queryIndex__model_activity_molecular_function_closure", curie=GOCAM.curie('model_activity_molecular_function_closure'),
                   model_uri=GOCAM.queryIndex__model_activity_molecular_function_closure, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_activity_molecular_function_rollup = Slot(uri=GOCAM.model_activity_molecular_function_rollup, name="queryIndex__model_activity_molecular_function_rollup", curie=GOCAM.curie('model_activity_molecular_function_rollup'),
                   model_uri=GOCAM.queryIndex__model_activity_molecular_function_rollup, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_activity_occurs_in_terms = Slot(uri=GOCAM.model_activity_occurs_in_terms, name="queryIndex__model_activity_occurs_in_terms", curie=GOCAM.curie('model_activity_occurs_in_terms'),
                   model_uri=GOCAM.queryIndex__model_activity_occurs_in_terms, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_activity_occurs_in_closure = Slot(uri=GOCAM.model_activity_occurs_in_closure, name="queryIndex__model_activity_occurs_in_closure", curie=GOCAM.curie('model_activity_occurs_in_closure'),
                   model_uri=GOCAM.queryIndex__model_activity_occurs_in_closure, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_activity_occurs_in_rollup = Slot(uri=GOCAM.model_activity_occurs_in_rollup, name="queryIndex__model_activity_occurs_in_rollup", curie=GOCAM.curie('model_activity_occurs_in_rollup'),
                   model_uri=GOCAM.queryIndex__model_activity_occurs_in_rollup, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_activity_enabled_by_terms = Slot(uri=GOCAM.model_activity_enabled_by_terms, name="queryIndex__model_activity_enabled_by_terms", curie=GOCAM.curie('model_activity_enabled_by_terms'),
                   model_uri=GOCAM.queryIndex__model_activity_enabled_by_terms, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_activity_enabled_by_closure = Slot(uri=GOCAM.model_activity_enabled_by_closure, name="queryIndex__model_activity_enabled_by_closure", curie=GOCAM.curie('model_activity_enabled_by_closure'),
                   model_uri=GOCAM.queryIndex__model_activity_enabled_by_closure, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_activity_enabled_by_rollup = Slot(uri=GOCAM.model_activity_enabled_by_rollup, name="queryIndex__model_activity_enabled_by_rollup", curie=GOCAM.curie('model_activity_enabled_by_rollup'),
                   model_uri=GOCAM.queryIndex__model_activity_enabled_by_rollup, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_activity_part_of_terms = Slot(uri=GOCAM.model_activity_part_of_terms, name="queryIndex__model_activity_part_of_terms", curie=GOCAM.curie('model_activity_part_of_terms'),
                   model_uri=GOCAM.queryIndex__model_activity_part_of_terms, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_activity_part_of_closure = Slot(uri=GOCAM.model_activity_part_of_closure, name="queryIndex__model_activity_part_of_closure", curie=GOCAM.curie('model_activity_part_of_closure'),
                   model_uri=GOCAM.queryIndex__model_activity_part_of_closure, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_activity_part_of_rollup = Slot(uri=GOCAM.model_activity_part_of_rollup, name="queryIndex__model_activity_part_of_rollup", curie=GOCAM.curie('model_activity_part_of_rollup'),
                   model_uri=GOCAM.queryIndex__model_activity_part_of_rollup, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_activity_has_input_terms = Slot(uri=GOCAM.model_activity_has_input_terms, name="queryIndex__model_activity_has_input_terms", curie=GOCAM.curie('model_activity_has_input_terms'),
                   model_uri=GOCAM.queryIndex__model_activity_has_input_terms, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_activity_has_input_closure = Slot(uri=GOCAM.model_activity_has_input_closure, name="queryIndex__model_activity_has_input_closure", curie=GOCAM.curie('model_activity_has_input_closure'),
                   model_uri=GOCAM.queryIndex__model_activity_has_input_closure, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_activity_has_input_rollup = Slot(uri=GOCAM.model_activity_has_input_rollup, name="queryIndex__model_activity_has_input_rollup", curie=GOCAM.curie('model_activity_has_input_rollup'),
                   model_uri=GOCAM.queryIndex__model_activity_has_input_rollup, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_taxon = Slot(uri=GOCAM.model_taxon, name="queryIndex__model_taxon", curie=GOCAM.curie('model_taxon'),
                   model_uri=GOCAM.queryIndex__model_taxon, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_taxon_closure = Slot(uri=GOCAM.model_taxon_closure, name="queryIndex__model_taxon_closure", curie=GOCAM.curie('model_taxon_closure'),
                   model_uri=GOCAM.queryIndex__model_taxon_closure, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__model_taxon_rollup = Slot(uri=GOCAM.model_taxon_rollup, name="queryIndex__model_taxon_rollup", curie=GOCAM.curie('model_taxon_rollup'),
                   model_uri=GOCAM.queryIndex__model_taxon_rollup, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__annoton_terms = Slot(uri=GOCAM.annoton_terms, name="queryIndex__annoton_terms", curie=GOCAM.curie('annoton_terms'),
                   model_uri=GOCAM.queryIndex__annoton_terms, domain=None, range=Optional[Union[dict[Union[str, TermObjectId], Union[dict, TermObject]], list[Union[dict, TermObject]]]])

slots.queryIndex__start_activities = Slot(uri=GOCAM.start_activities, name="queryIndex__start_activities", curie=GOCAM.curie('start_activities'),
                   model_uri=GOCAM.queryIndex__start_activities, domain=None, range=Optional[Union[Union[str, ActivityId], list[Union[str, ActivityId]]]])

slots.queryIndex__end_activities = Slot(uri=GOCAM.end_activities, name="queryIndex__end_activities", curie=GOCAM.curie('end_activities'),
                   model_uri=GOCAM.queryIndex__end_activities, domain=None, range=Optional[Union[Union[str, ActivityId], list[Union[str, ActivityId]]]])

slots.queryIndex__intermediate_activities = Slot(uri=GOCAM.intermediate_activities, name="queryIndex__intermediate_activities", curie=GOCAM.curie('intermediate_activities'),
                   model_uri=GOCAM.queryIndex__intermediate_activities, domain=None, range=Optional[Union[Union[str, ActivityId], list[Union[str, ActivityId]]]])

slots.queryIndex__singleton_activities = Slot(uri=GOCAM.singleton_activities, name="queryIndex__singleton_activities", curie=GOCAM.curie('singleton_activities'),
                   model_uri=GOCAM.queryIndex__singleton_activities, domain=None, range=Optional[Union[Union[str, ActivityId], list[Union[str, ActivityId]]]])

slots.queryIndex__number_of_start_activities = Slot(uri=GOCAM.number_of_start_activities, name="queryIndex__number_of_start_activities", curie=GOCAM.curie('number_of_start_activities'),
                   model_uri=GOCAM.queryIndex__number_of_start_activities, domain=None, range=Optional[int])

slots.queryIndex__number_of_end_activities = Slot(uri=GOCAM.number_of_end_activities, name="queryIndex__number_of_end_activities", curie=GOCAM.curie('number_of_end_activities'),
                   model_uri=GOCAM.queryIndex__number_of_end_activities, domain=None, range=Optional[int])

slots.queryIndex__number_of_intermediate_activities = Slot(uri=GOCAM.number_of_intermediate_activities, name="queryIndex__number_of_intermediate_activities", curie=GOCAM.curie('number_of_intermediate_activities'),
                   model_uri=GOCAM.queryIndex__number_of_intermediate_activities, domain=None, range=Optional[int])

slots.queryIndex__number_of_singleton_activities = Slot(uri=GOCAM.number_of_singleton_activities, name="queryIndex__number_of_singleton_activities", curie=GOCAM.curie('number_of_singleton_activities'),
                   model_uri=GOCAM.queryIndex__number_of_singleton_activities, domain=None, range=Optional[int])

slots.EnabledByGeneProductAssociation_term = Slot(uri=GOCAM.term, name="EnabledByGeneProductAssociation_term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.EnabledByGeneProductAssociation_term, domain=EnabledByGeneProductAssociation, range=Optional[Union[str, GeneProductTermObjectId]])

slots.EnabledByProteinComplexAssociation_term = Slot(uri=GOCAM.term, name="EnabledByProteinComplexAssociation_term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.EnabledByProteinComplexAssociation_term, domain=EnabledByProteinComplexAssociation, range=Optional[Union[str, ProteinComplexTermObjectId]])

slots.MolecularFunctionAssociation_term = Slot(uri=GOCAM.term, name="MolecularFunctionAssociation_term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.MolecularFunctionAssociation_term, domain=MolecularFunctionAssociation, range=Optional[Union[str, MolecularFunctionTermObjectId]])

slots.BiologicalProcessAssociation_term = Slot(uri=GOCAM.term, name="BiologicalProcessAssociation_term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.BiologicalProcessAssociation_term, domain=BiologicalProcessAssociation, range=Optional[Union[str, BiologicalProcessTermObjectId]])

slots.CellularAnatomicalEntityAssociation_term = Slot(uri=GOCAM.term, name="CellularAnatomicalEntityAssociation_term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.CellularAnatomicalEntityAssociation_term, domain=CellularAnatomicalEntityAssociation, range=Optional[Union[str, CellularAnatomicalEntityTermObjectId]])

slots.CellTypeAssociation_term = Slot(uri=GOCAM.term, name="CellTypeAssociation_term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.CellTypeAssociation_term, domain=CellTypeAssociation, range=Optional[Union[str, CellTypeTermObjectId]])

slots.GrossAnatomyAssociation_term = Slot(uri=GOCAM.term, name="GrossAnatomyAssociation_term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.GrossAnatomyAssociation_term, domain=GrossAnatomyAssociation, range=Optional[Union[str, GrossAnatomicalStructureTermObjectId]])

slots.MoleculeAssociation_term = Slot(uri=GOCAM.term, name="MoleculeAssociation_term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.MoleculeAssociation_term, domain=MoleculeAssociation, range=Optional[Union[str, MoleculeTermObjectId]])