# Auto generated from gocam.yaml by pythongen.py version: 0.0.1
# Generation date: 2024-08-28T17:21:26
# Schema: gocam
#
# id: https://w3id.org/gocam
# description: GO CAM LinkML schema (experimental)
#
#   The central class in this datamodel is a [Model](Model.md). A model consists of a set of
#   [Activity](Activity.md) objects.
# license: https://creativecommons.org/publicdomain/zero/1.0/

import dataclasses
import re
from jsonasobj2 import JsonObj, as_dict
from typing import Optional, List, Union, Dict, ClassVar, Any
from dataclasses import dataclass
from datetime import date, datetime
from linkml_runtime.linkml_model.meta import EnumDefinition, PermissibleValue, PvFormulaOptions

from linkml_runtime.utils.slot import Slot
from linkml_runtime.utils.metamodelcore import empty_list, empty_dict, bnode
from linkml_runtime.utils.yamlutils import YAMLRoot, extended_str, extended_float, extended_int
from linkml_runtime.utils.dataclass_extensions_376 import dataclasses_init_fn_with_kwargs
from linkml_runtime.utils.formatutils import camelcase, underscore, sfx
from linkml_runtime.utils.enumerations import EnumDefinitionImpl
from rdflib import Namespace, URIRef
from linkml_runtime.utils.curienamespace import CurieNamespace
from linkml_runtime.linkml_model.types import Boolean, String, Uriorcurie
from linkml_runtime.utils.metamodelcore import Bool, URIorCURIE

metamodel_version = "1.7.0"
version = None

# Overwrite dataclasses _init_fn to add **kwargs in __init__
dataclasses._init_fn = dataclasses_init_fn_with_kwargs

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
PMID = CurieNamespace('PMID', 'http://identifiers.org/pmid/')
PO = CurieNamespace('PO', 'http://example.org/UNKNOWN/PO/')
RO = CurieNamespace('RO', 'http://purl.obolibrary.org/obo/RO_')
UBERON = CurieNamespace('UBERON', 'http://example.org/UNKNOWN/UBERON/')
UNIPROTKB = CurieNamespace('UniProtKB', 'http://identifiers.org/uniprot/')
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


@dataclass
class Model(YAMLRoot):
    """
    A model of a biological program consisting of a set of causally connected activities
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["Model"]
    class_class_curie: ClassVar[str] = "gocam:Model"
    class_name: ClassVar[str] = "Model"
    class_model_uri: ClassVar[URIRef] = GOCAM.Model

    id: Union[str, ModelId] = None
    title: Optional[str] = None
    taxon: Optional[Union[str, TaxonTermObjectId]] = None
    status: Optional[Union[str, "ModelStateEnum"]] = None
    comments: Optional[Union[str, List[str]]] = empty_list()
    activities: Optional[Union[Dict[Union[str, ActivityId], Union[dict, "Activity"]], List[Union[dict, "Activity"]]]] = empty_dict()
    objects: Optional[Union[Dict[Union[str, ObjectId], Union[dict, "Object"]], List[Union[dict, "Object"]]]] = empty_dict()
    provenances: Optional[Union[Union[dict, "ProvenanceInfo"], List[Union[dict, "ProvenanceInfo"]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, ModelId):
            self.id = ModelId(self.id)

        if self.title is not None and not isinstance(self.title, str):
            self.title = str(self.title)

        if self.taxon is not None and not isinstance(self.taxon, TaxonTermObjectId):
            self.taxon = TaxonTermObjectId(self.taxon)

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

        super().__post_init__(**kwargs)


@dataclass
class Activity(YAMLRoot):
    """
    An individual activity in a causal model, representing the individual molecular activity of a single gene product
    or complex in the context of a particular model
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["Activity"]
    class_class_curie: ClassVar[str] = "gocam:Activity"
    class_name: ClassVar[str] = "Activity"
    class_model_uri: ClassVar[URIRef] = GOCAM.Activity

    id: Union[str, ActivityId] = None
    enabled_by: Optional[Union[dict, "EnabledByAssociation"]] = None
    molecular_function: Optional[Union[dict, "MolecularFunctionAssociation"]] = None
    occurs_in: Optional[Union[dict, "CellularAnatomicalEntityAssociation"]] = None
    part_of: Optional[Union[dict, "BiologicalProcessAssociation"]] = None
    has_direct_input: Optional[Union[dict, "MoleculeAssociation"]] = None
    causal_associations: Optional[Union[Union[dict, "CausalAssociation"], List[Union[dict, "CausalAssociation"]]]] = empty_list()
    provenances: Optional[Union[Union[dict, "ProvenanceInfo"], List[Union[dict, "ProvenanceInfo"]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
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

        if self.has_direct_input is not None and not isinstance(self.has_direct_input, MoleculeAssociation):
            self.has_direct_input = MoleculeAssociation(**as_dict(self.has_direct_input))

        if not isinstance(self.causal_associations, list):
            self.causal_associations = [self.causal_associations] if self.causal_associations is not None else []
        self.causal_associations = [v if isinstance(v, CausalAssociation) else CausalAssociation(**as_dict(v)) for v in self.causal_associations]

        if not isinstance(self.provenances, list):
            self.provenances = [self.provenances] if self.provenances is not None else []
        self.provenances = [v if isinstance(v, ProvenanceInfo) else ProvenanceInfo(**as_dict(v)) for v in self.provenances]

        super().__post_init__(**kwargs)


@dataclass
class EvidenceItem(YAMLRoot):
    """
    An individual piece of evidence that is associated with an assertion in a model
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["EvidenceItem"]
    class_class_curie: ClassVar[str] = "gocam:EvidenceItem"
    class_name: ClassVar[str] = "EvidenceItem"
    class_model_uri: ClassVar[URIRef] = GOCAM.EvidenceItem

    term: Optional[Union[str, EvidenceTermObjectId]] = None
    reference: Optional[Union[str, PublicationObjectId]] = None
    with_objects: Optional[Union[Union[str, ObjectId], List[Union[str, ObjectId]]]] = empty_list()
    provenances: Optional[Union[Union[dict, "ProvenanceInfo"], List[Union[dict, "ProvenanceInfo"]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
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


@dataclass
class Association(YAMLRoot):
    """
    An abstract grouping for different kinds of evidence-associated provenance
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["Association"]
    class_class_curie: ClassVar[str] = "gocam:Association"
    class_name: ClassVar[str] = "Association"
    class_model_uri: ClassVar[URIRef] = GOCAM.Association

    type: Optional[str] = None
    evidence: Optional[Union[Union[dict, EvidenceItem], List[Union[dict, EvidenceItem]]]] = empty_list()
    provenances: Optional[Union[Union[dict, "ProvenanceInfo"], List[Union[dict, "ProvenanceInfo"]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
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



@dataclass
class EnabledByAssociation(Association):
    """
    An association between an activity and the gene product or complex that carries it out
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["EnabledByAssociation"]
    class_class_curie: ClassVar[str] = "gocam:EnabledByAssociation"
    class_name: ClassVar[str] = "EnabledByAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.EnabledByAssociation

    term: Optional[Union[str, InformationBiomacromoleculeTermObjectId]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.term is not None and not isinstance(self.term, InformationBiomacromoleculeTermObjectId):
            self.term = InformationBiomacromoleculeTermObjectId(self.term)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass
class EnabledByGeneProductAssociation(EnabledByAssociation):
    """
    An association between an activity and a gene product
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["EnabledByGeneProductAssociation"]
    class_class_curie: ClassVar[str] = "gocam:EnabledByGeneProductAssociation"
    class_name: ClassVar[str] = "EnabledByGeneProductAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.EnabledByGeneProductAssociation

    term: Optional[Union[str, GeneProductTermObjectId]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.term is not None and not isinstance(self.term, GeneProductTermObjectId):
            self.term = GeneProductTermObjectId(self.term)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass
class EnabledByProteinComplexAssociation(EnabledByAssociation):
    """
    An association between an activity and a protein complex
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["EnabledByProteinComplexAssociation"]
    class_class_curie: ClassVar[str] = "gocam:EnabledByProteinComplexAssociation"
    class_name: ClassVar[str] = "EnabledByProteinComplexAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.EnabledByProteinComplexAssociation

    members: Optional[Union[Union[str, GeneProductTermObjectId], List[Union[str, GeneProductTermObjectId]]]] = empty_list()
    term: Optional[Union[str, ProteinComplexTermObjectId]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if not isinstance(self.members, list):
            self.members = [self.members] if self.members is not None else []
        self.members = [v if isinstance(v, GeneProductTermObjectId) else GeneProductTermObjectId(v) for v in self.members]

        if self.term is not None and not isinstance(self.term, ProteinComplexTermObjectId):
            self.term = ProteinComplexTermObjectId(self.term)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass
class CausalAssociation(Association):
    """
    A causal association between two activities
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["CausalAssociation"]
    class_class_curie: ClassVar[str] = "gocam:CausalAssociation"
    class_name: ClassVar[str] = "CausalAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.CausalAssociation

    predicate: Optional[Union[str, PredicateTermObjectId]] = None
    downstream_activity: Optional[Union[str, ActivityId]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.predicate is not None and not isinstance(self.predicate, PredicateTermObjectId):
            self.predicate = PredicateTermObjectId(self.predicate)

        if self.downstream_activity is not None and not isinstance(self.downstream_activity, ActivityId):
            self.downstream_activity = ActivityId(self.downstream_activity)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass
class TermAssociation(Association):
    """
    An association between an activity and a term, potentially with extensions
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["TermAssociation"]
    class_class_curie: ClassVar[str] = "gocam:TermAssociation"
    class_name: ClassVar[str] = "TermAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.TermAssociation

    term: Optional[Union[str, TermObjectId]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.term is not None and not isinstance(self.term, TermObjectId):
            self.term = TermObjectId(self.term)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass
class MolecularFunctionAssociation(TermAssociation):
    """
    An association between an activity and a molecular function term
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["MolecularFunctionAssociation"]
    class_class_curie: ClassVar[str] = "gocam:MolecularFunctionAssociation"
    class_name: ClassVar[str] = "MolecularFunctionAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.MolecularFunctionAssociation

    term: Optional[Union[str, MolecularFunctionTermObjectId]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.term is not None and not isinstance(self.term, MolecularFunctionTermObjectId):
            self.term = MolecularFunctionTermObjectId(self.term)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass
class BiologicalProcessAssociation(TermAssociation):
    """
    An association between an activity and a biological process term
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["BiologicalProcessAssociation"]
    class_class_curie: ClassVar[str] = "gocam:BiologicalProcessAssociation"
    class_name: ClassVar[str] = "BiologicalProcessAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.BiologicalProcessAssociation

    term: Optional[Union[str, BiologicalProcessTermObjectId]] = None
    happens_during: Optional[Union[str, PhaseTermObjectId]] = None
    part_of: Optional[Union[str, BiologicalProcessTermObjectId]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.term is not None and not isinstance(self.term, BiologicalProcessTermObjectId):
            self.term = BiologicalProcessTermObjectId(self.term)

        if self.happens_during is not None and not isinstance(self.happens_during, PhaseTermObjectId):
            self.happens_during = PhaseTermObjectId(self.happens_during)

        if self.part_of is not None and not isinstance(self.part_of, BiologicalProcessTermObjectId):
            self.part_of = BiologicalProcessTermObjectId(self.part_of)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass
class CellularAnatomicalEntityAssociation(TermAssociation):
    """
    An association between an activity and a cellular anatomical entity term
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["CellularAnatomicalEntityAssociation"]
    class_class_curie: ClassVar[str] = "gocam:CellularAnatomicalEntityAssociation"
    class_name: ClassVar[str] = "CellularAnatomicalEntityAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.CellularAnatomicalEntityAssociation

    term: Optional[Union[str, CellularAnatomicalEntityTermObjectId]] = None
    part_of: Optional[Union[dict, "CellTypeAssociation"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.term is not None and not isinstance(self.term, CellularAnatomicalEntityTermObjectId):
            self.term = CellularAnatomicalEntityTermObjectId(self.term)

        if self.part_of is not None and not isinstance(self.part_of, CellTypeAssociation):
            self.part_of = CellTypeAssociation(**as_dict(self.part_of))

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass
class CellTypeAssociation(TermAssociation):
    """
    An association between an activity and a cell type term
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["CellTypeAssociation"]
    class_class_curie: ClassVar[str] = "gocam:CellTypeAssociation"
    class_name: ClassVar[str] = "CellTypeAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.CellTypeAssociation

    term: Optional[Union[str, CellTypeTermObjectId]] = None
    part_of: Optional[Union[dict, "GrossAnatomyAssociation"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.term is not None and not isinstance(self.term, CellTypeTermObjectId):
            self.term = CellTypeTermObjectId(self.term)

        if self.part_of is not None and not isinstance(self.part_of, GrossAnatomyAssociation):
            self.part_of = GrossAnatomyAssociation(**as_dict(self.part_of))

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass
class GrossAnatomyAssociation(TermAssociation):
    """
    An association between an activity and a gross anatomical structure term
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["GrossAnatomyAssociation"]
    class_class_curie: ClassVar[str] = "gocam:GrossAnatomyAssociation"
    class_name: ClassVar[str] = "GrossAnatomyAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.GrossAnatomyAssociation

    term: Optional[Union[str, GrossAnatomicalStructureTermObjectId]] = None
    part_of: Optional[Union[dict, "GrossAnatomyAssociation"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.term is not None and not isinstance(self.term, GrossAnatomicalStructureTermObjectId):
            self.term = GrossAnatomicalStructureTermObjectId(self.term)

        if self.part_of is not None and not isinstance(self.part_of, GrossAnatomyAssociation):
            self.part_of = GrossAnatomyAssociation(**as_dict(self.part_of))

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass
class MoleculeAssociation(TermAssociation):
    """
    An association between an activity and a molecule term
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["MoleculeAssociation"]
    class_class_curie: ClassVar[str] = "gocam:MoleculeAssociation"
    class_name: ClassVar[str] = "MoleculeAssociation"
    class_model_uri: ClassVar[URIRef] = GOCAM.MoleculeAssociation

    term: Optional[Union[str, MoleculeTermObjectId]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.term is not None and not isinstance(self.term, MoleculeTermObjectId):
            self.term = MoleculeTermObjectId(self.term)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_name)


@dataclass
class Object(YAMLRoot):
    """
    An abstract class for all identified objects in a model
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["Object"]
    class_class_curie: ClassVar[str] = "gocam:Object"
    class_name: ClassVar[str] = "Object"
    class_model_uri: ClassVar[URIRef] = GOCAM.Object

    id: Union[str, ObjectId] = None
    label: Optional[str] = None
    type: Optional[Union[str, URIorCURIE]] = None
    obsolete: Optional[Union[bool, Bool]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
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



@dataclass
class TermObject(Object):
    """
    An abstract class for all ontology term objects
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["TermObject"]
    class_class_curie: ClassVar[str] = "gocam:TermObject"
    class_name: ClassVar[str] = "TermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.TermObject

    id: Union[str, TermObjectId] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass
class PublicationObject(Object):
    """
    An object that represents a publication or other kind of reference
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["PublicationObject"]
    class_class_curie: ClassVar[str] = "gocam:PublicationObject"
    class_name: ClassVar[str] = "PublicationObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.PublicationObject

    id: Union[str, PublicationObjectId] = None
    abstract_text: Optional[str] = None
    full_text: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
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


@dataclass
class EvidenceTermObject(TermObject):
    """
    A term object that represents an evidence term from ECO
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["EvidenceTermObject"]
    class_class_curie: ClassVar[str] = "gocam:EvidenceTermObject"
    class_name: ClassVar[str] = "EvidenceTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.EvidenceTermObject

    id: Union[str, EvidenceTermObjectId] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, EvidenceTermObjectId):
            self.id = EvidenceTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass
class MolecularFunctionTermObject(TermObject):
    """
    A term object that represents a molecular function term from GO
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["MolecularFunctionTermObject"]
    class_class_curie: ClassVar[str] = "gocam:MolecularFunctionTermObject"
    class_name: ClassVar[str] = "MolecularFunctionTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.MolecularFunctionTermObject

    id: Union[str, MolecularFunctionTermObjectId] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, MolecularFunctionTermObjectId):
            self.id = MolecularFunctionTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass
class BiologicalProcessTermObject(TermObject):
    """
    A term object that represents a biological process term from GO
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["BiologicalProcessTermObject"]
    class_class_curie: ClassVar[str] = "gocam:BiologicalProcessTermObject"
    class_name: ClassVar[str] = "BiologicalProcessTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.BiologicalProcessTermObject

    id: Union[str, BiologicalProcessTermObjectId] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, BiologicalProcessTermObjectId):
            self.id = BiologicalProcessTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass
class CellularAnatomicalEntityTermObject(TermObject):
    """
    A term object that represents a cellular anatomical entity term from GO
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["CellularAnatomicalEntityTermObject"]
    class_class_curie: ClassVar[str] = "gocam:CellularAnatomicalEntityTermObject"
    class_name: ClassVar[str] = "CellularAnatomicalEntityTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.CellularAnatomicalEntityTermObject

    id: Union[str, CellularAnatomicalEntityTermObjectId] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, CellularAnatomicalEntityTermObjectId):
            self.id = CellularAnatomicalEntityTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass
class MoleculeTermObject(TermObject):
    """
    A term object that represents a molecule term from CHEBI or UniProtKB
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["MoleculeTermObject"]
    class_class_curie: ClassVar[str] = "gocam:MoleculeTermObject"
    class_name: ClassVar[str] = "MoleculeTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.MoleculeTermObject

    id: Union[str, MoleculeTermObjectId] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, MoleculeTermObjectId):
            self.id = MoleculeTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass
class CellTypeTermObject(TermObject):
    """
    A term object that represents a cell type term from CL
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["CellTypeTermObject"]
    class_class_curie: ClassVar[str] = "gocam:CellTypeTermObject"
    class_name: ClassVar[str] = "CellTypeTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.CellTypeTermObject

    id: Union[str, CellTypeTermObjectId] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, CellTypeTermObjectId):
            self.id = CellTypeTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass
class GrossAnatomicalStructureTermObject(TermObject):
    """
    A term object that represents a gross anatomical structure term from UBERON
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["GrossAnatomicalStructureTermObject"]
    class_class_curie: ClassVar[str] = "gocam:GrossAnatomicalStructureTermObject"
    class_name: ClassVar[str] = "GrossAnatomicalStructureTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.GrossAnatomicalStructureTermObject

    id: Union[str, GrossAnatomicalStructureTermObjectId] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, GrossAnatomicalStructureTermObjectId):
            self.id = GrossAnatomicalStructureTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass
class PhaseTermObject(TermObject):
    """
    A term object that represents a phase term from GO or UBERON
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["PhaseTermObject"]
    class_class_curie: ClassVar[str] = "gocam:PhaseTermObject"
    class_name: ClassVar[str] = "PhaseTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.PhaseTermObject

    id: Union[str, PhaseTermObjectId] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, PhaseTermObjectId):
            self.id = PhaseTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass
class InformationBiomacromoleculeTermObject(TermObject):
    """
    An abstract class for all information biomacromolecule term objects
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["InformationBiomacromoleculeTermObject"]
    class_class_curie: ClassVar[str] = "gocam:InformationBiomacromoleculeTermObject"
    class_name: ClassVar[str] = "InformationBiomacromoleculeTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.InformationBiomacromoleculeTermObject

    id: Union[str, InformationBiomacromoleculeTermObjectId] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass
class GeneProductTermObject(InformationBiomacromoleculeTermObject):
    """
    A term object that represents a gene product term from GO or UniProtKB
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["GeneProductTermObject"]
    class_class_curie: ClassVar[str] = "gocam:GeneProductTermObject"
    class_name: ClassVar[str] = "GeneProductTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.GeneProductTermObject

    id: Union[str, GeneProductTermObjectId] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, GeneProductTermObjectId):
            self.id = GeneProductTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass
class ProteinComplexTermObject(InformationBiomacromoleculeTermObject):
    """
    A term object that represents a protein complex term from GO
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["ProteinComplexTermObject"]
    class_class_curie: ClassVar[str] = "gocam:ProteinComplexTermObject"
    class_name: ClassVar[str] = "ProteinComplexTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.ProteinComplexTermObject

    id: Union[str, ProteinComplexTermObjectId] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, ProteinComplexTermObjectId):
            self.id = ProteinComplexTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass
class TaxonTermObject(TermObject):
    """
    A term object that represents a taxon term from NCBITaxon
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["TaxonTermObject"]
    class_class_curie: ClassVar[str] = "gocam:TaxonTermObject"
    class_name: ClassVar[str] = "TaxonTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.TaxonTermObject

    id: Union[str, TaxonTermObjectId] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, TaxonTermObjectId):
            self.id = TaxonTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass
class PredicateTermObject(TermObject):
    """
    A term object that represents a taxon term from NCBITaxon
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["PredicateTermObject"]
    class_class_curie: ClassVar[str] = "gocam:PredicateTermObject"
    class_name: ClassVar[str] = "PredicateTermObject"
    class_model_uri: ClassVar[URIRef] = GOCAM.PredicateTermObject

    id: Union[str, PredicateTermObjectId] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, PredicateTermObjectId):
            self.id = PredicateTermObjectId(self.id)

        super().__post_init__(**kwargs)
        self.unknown_type = str(self.class_class_curie)


@dataclass
class ProvenanceInfo(YAMLRoot):
    """
    Provenance information for an object
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GOCAM["ProvenanceInfo"]
    class_class_curie: ClassVar[str] = "gocam:ProvenanceInfo"
    class_name: ClassVar[str] = "ProvenanceInfo"
    class_model_uri: ClassVar[URIRef] = GOCAM.ProvenanceInfo

    contributor: Optional[str] = None
    created: Optional[str] = None
    date: Optional[str] = None
    provided_by: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.contributor is not None and not isinstance(self.contributor, str):
            self.contributor = str(self.contributor)

        if self.created is not None and not isinstance(self.created, str):
            self.created = str(self.created)

        if self.date is not None and not isinstance(self.date, str):
            self.date = str(self.date)

        if self.provided_by is not None and not isinstance(self.provided_by, str):
            self.provided_by = str(self.provided_by)

        super().__post_init__(**kwargs)


# Enumerations
class ModelStateEnum(EnumDefinitionImpl):
    """
    Status of a model
    """
    production = PermissibleValue(text="production")
    development = PermissibleValue(text="development")

    _defn = EnumDefinition(
        name="ModelStateEnum",
        description="Status of a model",
    )

class InformationBiomacromoleculeCategory(EnumDefinitionImpl):

    GeneOrReferenceProtein = PermissibleValue(
        text="GeneOrReferenceProtein",
        meaning=GOCAM["biolink.GeneOrGeneProduct"])
    ProteinIsoform = PermissibleValue(text="ProteinIsoform")
    MacromolecularComplex = PermissibleValue(text="MacromolecularComplex")
    Unknown = PermissibleValue(text="Unknown")

    _defn = EnumDefinition(
        name="InformationBiomacromoleculeCategory",
    )

class CausalPredicateEnum(EnumDefinitionImpl):

    regulates = PermissibleValue(
        text="regulates",
        meaning=RO["0002211"])

    _defn = EnumDefinition(
        name="CausalPredicateEnum",
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

# Slots
class slots:
    pass

slots.model__id = Slot(uri=GOCAM.id, name="model__id", curie=GOCAM.curie('id'),
                   model_uri=GOCAM.model__id, domain=None, range=URIRef)

slots.model__title = Slot(uri=DCT.title, name="model__title", curie=DCT.curie('title'),
                   model_uri=GOCAM.model__title, domain=None, range=Optional[str])

slots.model__taxon = Slot(uri=GOCAM.taxon, name="model__taxon", curie=GOCAM.curie('taxon'),
                   model_uri=GOCAM.model__taxon, domain=None, range=Optional[Union[str, TaxonTermObjectId]])

slots.model__status = Slot(uri=PAV.status, name="model__status", curie=PAV.curie('status'),
                   model_uri=GOCAM.model__status, domain=None, range=Optional[Union[str, "ModelStateEnum"]])

slots.model__comments = Slot(uri=RDFS.comment, name="model__comments", curie=RDFS.curie('comment'),
                   model_uri=GOCAM.model__comments, domain=None, range=Optional[Union[str, List[str]]])

slots.model__activities = Slot(uri=GOCAM.activities, name="model__activities", curie=GOCAM.curie('activities'),
                   model_uri=GOCAM.model__activities, domain=None, range=Optional[Union[Dict[Union[str, ActivityId], Union[dict, Activity]], List[Union[dict, Activity]]]])

slots.model__objects = Slot(uri=GOCAM.objects, name="model__objects", curie=GOCAM.curie('objects'),
                   model_uri=GOCAM.model__objects, domain=None, range=Optional[Union[Dict[Union[str, ObjectId], Union[dict, Object]], List[Union[dict, Object]]]])

slots.model__provenances = Slot(uri=GOCAM.provenances, name="model__provenances", curie=GOCAM.curie('provenances'),
                   model_uri=GOCAM.model__provenances, domain=None, range=Optional[Union[Union[dict, ProvenanceInfo], List[Union[dict, ProvenanceInfo]]]])

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

slots.activity__has_direct_input = Slot(uri=GOCAM.has_direct_input, name="activity__has_direct_input", curie=GOCAM.curie('has_direct_input'),
                   model_uri=GOCAM.activity__has_direct_input, domain=None, range=Optional[Union[dict, MoleculeAssociation]])

slots.activity__causal_associations = Slot(uri=GOCAM.causal_associations, name="activity__causal_associations", curie=GOCAM.curie('causal_associations'),
                   model_uri=GOCAM.activity__causal_associations, domain=None, range=Optional[Union[Union[dict, CausalAssociation], List[Union[dict, CausalAssociation]]]])

slots.activity__provenances = Slot(uri=GOCAM.provenances, name="activity__provenances", curie=GOCAM.curie('provenances'),
                   model_uri=GOCAM.activity__provenances, domain=None, range=Optional[Union[Union[dict, ProvenanceInfo], List[Union[dict, ProvenanceInfo]]]])

slots.evidenceItem__term = Slot(uri=GOCAM.term, name="evidenceItem__term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.evidenceItem__term, domain=None, range=Optional[Union[str, EvidenceTermObjectId]])

slots.evidenceItem__reference = Slot(uri=GOCAM.reference, name="evidenceItem__reference", curie=GOCAM.curie('reference'),
                   model_uri=GOCAM.evidenceItem__reference, domain=None, range=Optional[Union[str, PublicationObjectId]])

slots.evidenceItem__with_objects = Slot(uri=GOCAM.with_objects, name="evidenceItem__with_objects", curie=GOCAM.curie('with_objects'),
                   model_uri=GOCAM.evidenceItem__with_objects, domain=None, range=Optional[Union[Union[str, ObjectId], List[Union[str, ObjectId]]]])

slots.evidenceItem__provenances = Slot(uri=GOCAM.provenances, name="evidenceItem__provenances", curie=GOCAM.curie('provenances'),
                   model_uri=GOCAM.evidenceItem__provenances, domain=None, range=Optional[Union[Union[dict, ProvenanceInfo], List[Union[dict, ProvenanceInfo]]]])

slots.association__type = Slot(uri=GOCAM.type, name="association__type", curie=GOCAM.curie('type'),
                   model_uri=GOCAM.association__type, domain=None, range=Optional[str])

slots.association__evidence = Slot(uri=GOCAM.evidence, name="association__evidence", curie=GOCAM.curie('evidence'),
                   model_uri=GOCAM.association__evidence, domain=None, range=Optional[Union[Union[dict, EvidenceItem], List[Union[dict, EvidenceItem]]]])

slots.association__provenances = Slot(uri=GOCAM.provenances, name="association__provenances", curie=GOCAM.curie('provenances'),
                   model_uri=GOCAM.association__provenances, domain=None, range=Optional[Union[Union[dict, ProvenanceInfo], List[Union[dict, ProvenanceInfo]]]])

slots.enabledByAssociation__term = Slot(uri=GOCAM.term, name="enabledByAssociation__term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.enabledByAssociation__term, domain=None, range=Optional[Union[str, InformationBiomacromoleculeTermObjectId]])

slots.enabledByProteinComplexAssociation__members = Slot(uri=GOCAM.members, name="enabledByProteinComplexAssociation__members", curie=GOCAM.curie('members'),
                   model_uri=GOCAM.enabledByProteinComplexAssociation__members, domain=None, range=Optional[Union[Union[str, GeneProductTermObjectId], List[Union[str, GeneProductTermObjectId]]]])

slots.causalAssociation__predicate = Slot(uri=GOCAM.predicate, name="causalAssociation__predicate", curie=GOCAM.curie('predicate'),
                   model_uri=GOCAM.causalAssociation__predicate, domain=None, range=Optional[Union[str, PredicateTermObjectId]])

slots.causalAssociation__downstream_activity = Slot(uri=GOCAM.downstream_activity, name="causalAssociation__downstream_activity", curie=GOCAM.curie('downstream_activity'),
                   model_uri=GOCAM.causalAssociation__downstream_activity, domain=None, range=Optional[Union[str, ActivityId]])

slots.termAssociation__term = Slot(uri=GOCAM.term, name="termAssociation__term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.termAssociation__term, domain=None, range=Optional[Union[str, TermObjectId]])

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
                   model_uri=GOCAM.provenanceInfo__contributor, domain=None, range=Optional[str])

slots.provenanceInfo__created = Slot(uri=DCT.created, name="provenanceInfo__created", curie=DCT.curie('created'),
                   model_uri=GOCAM.provenanceInfo__created, domain=None, range=Optional[str])

slots.provenanceInfo__date = Slot(uri=DCT.date, name="provenanceInfo__date", curie=DCT.curie('date'),
                   model_uri=GOCAM.provenanceInfo__date, domain=None, range=Optional[str])

slots.provenanceInfo__provided_by = Slot(uri=PAV.providedBy, name="provenanceInfo__provided_by", curie=PAV.curie('providedBy'),
                   model_uri=GOCAM.provenanceInfo__provided_by, domain=None, range=Optional[str])

slots.happens_during = Slot(uri=GOCAM.happens_during, name="happens_during", curie=GOCAM.curie('happens_during'),
                   model_uri=GOCAM.happens_during, domain=None, range=Optional[Union[str, PhaseTermObjectId]])

slots.part_of = Slot(uri=GOCAM.part_of, name="part_of", curie=GOCAM.curie('part_of'),
                   model_uri=GOCAM.part_of, domain=None, range=Optional[Union[str, BiologicalProcessTermObjectId]])

slots.EnabledByGeneProductAssociation_term = Slot(uri=GOCAM.term, name="EnabledByGeneProductAssociation_term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.EnabledByGeneProductAssociation_term, domain=EnabledByGeneProductAssociation, range=Optional[Union[str, GeneProductTermObjectId]])

slots.EnabledByProteinComplexAssociation_term = Slot(uri=GOCAM.term, name="EnabledByProteinComplexAssociation_term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.EnabledByProteinComplexAssociation_term, domain=EnabledByProteinComplexAssociation, range=Optional[Union[str, ProteinComplexTermObjectId]])

slots.MolecularFunctionAssociation_term = Slot(uri=GOCAM.term, name="MolecularFunctionAssociation_term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.MolecularFunctionAssociation_term, domain=MolecularFunctionAssociation, range=Optional[Union[str, MolecularFunctionTermObjectId]])

slots.BiologicalProcessAssociation_term = Slot(uri=GOCAM.term, name="BiologicalProcessAssociation_term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.BiologicalProcessAssociation_term, domain=BiologicalProcessAssociation, range=Optional[Union[str, BiologicalProcessTermObjectId]])

slots.BiologicalProcessAssociation_happens_during = Slot(uri=GOCAM.happens_during, name="BiologicalProcessAssociation_happens_during", curie=GOCAM.curie('happens_during'),
                   model_uri=GOCAM.BiologicalProcessAssociation_happens_during, domain=BiologicalProcessAssociation, range=Optional[Union[str, PhaseTermObjectId]])

slots.BiologicalProcessAssociation_part_of = Slot(uri=GOCAM.part_of, name="BiologicalProcessAssociation_part_of", curie=GOCAM.curie('part_of'),
                   model_uri=GOCAM.BiologicalProcessAssociation_part_of, domain=BiologicalProcessAssociation, range=Optional[Union[str, BiologicalProcessTermObjectId]])

slots.CellularAnatomicalEntityAssociation_term = Slot(uri=GOCAM.term, name="CellularAnatomicalEntityAssociation_term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.CellularAnatomicalEntityAssociation_term, domain=CellularAnatomicalEntityAssociation, range=Optional[Union[str, CellularAnatomicalEntityTermObjectId]])

slots.CellularAnatomicalEntityAssociation_part_of = Slot(uri=GOCAM.part_of, name="CellularAnatomicalEntityAssociation_part_of", curie=GOCAM.curie('part_of'),
                   model_uri=GOCAM.CellularAnatomicalEntityAssociation_part_of, domain=CellularAnatomicalEntityAssociation, range=Optional[Union[dict, "CellTypeAssociation"]])

slots.CellTypeAssociation_term = Slot(uri=GOCAM.term, name="CellTypeAssociation_term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.CellTypeAssociation_term, domain=CellTypeAssociation, range=Optional[Union[str, CellTypeTermObjectId]])

slots.CellTypeAssociation_part_of = Slot(uri=GOCAM.part_of, name="CellTypeAssociation_part_of", curie=GOCAM.curie('part_of'),
                   model_uri=GOCAM.CellTypeAssociation_part_of, domain=CellTypeAssociation, range=Optional[Union[dict, "GrossAnatomyAssociation"]])

slots.GrossAnatomyAssociation_term = Slot(uri=GOCAM.term, name="GrossAnatomyAssociation_term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.GrossAnatomyAssociation_term, domain=GrossAnatomyAssociation, range=Optional[Union[str, GrossAnatomicalStructureTermObjectId]])

slots.GrossAnatomyAssociation_part_of = Slot(uri=GOCAM.part_of, name="GrossAnatomyAssociation_part_of", curie=GOCAM.curie('part_of'),
                   model_uri=GOCAM.GrossAnatomyAssociation_part_of, domain=GrossAnatomyAssociation, range=Optional[Union[dict, "GrossAnatomyAssociation"]])

slots.MoleculeAssociation_term = Slot(uri=GOCAM.term, name="MoleculeAssociation_term", curie=GOCAM.curie('term'),
                   model_uri=GOCAM.MoleculeAssociation_term, domain=MoleculeAssociation, range=Optional[Union[str, MoleculeTermObjectId]])