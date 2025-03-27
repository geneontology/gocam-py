from __future__ import annotations 

import re
import sys
from datetime import (
    date,
    datetime,
    time
)
from decimal import Decimal 
from enum import Enum 
from typing import (
    Any,
    ClassVar,
    Dict,
    List,
    Literal,
    Optional,
    Union
)

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    RootModel,
    field_validator
)


metamodel_version = "None"
version = "None"


class ConfiguredBaseModel(BaseModel):
    model_config = ConfigDict(
        validate_assignment = True,
        validate_default = True,
        extra = "forbid",
        arbitrary_types_allowed = True,
        use_enum_values = True,
        strict = False,
    )
    pass




class LinkMLMeta(RootModel):
    root: Dict[str, Any] = {}
    model_config = ConfigDict(frozen=True)

    def __getattr__(self, key:str):
        return getattr(self.root, key)

    def __getitem__(self, key:str):
        return self.root[key]

    def __setitem__(self, key:str, value):
        self.root[key] = value

    def __contains__(self, key:str) -> bool:
        return key in self.root


linkml_meta = LinkMLMeta({'default_prefix': 'gocam',
     'default_range': 'string',
     'description': 'GO CAM LinkML schema (experimental)\n'
                    '\n'
                    'The central class in this datamodel is a [Model](Model.md). A '
                    'model consists of a set of\n'
                    '[Activity](Activity.md) objects.',
     'id': 'https://w3id.org/gocam',
     'imports': ['linkml:types'],
     'name': 'gocam',
     'prefixes': {'BFO': {'prefix_prefix': 'BFO',
                          'prefix_reference': 'http://purl.obolibrary.org/obo/BFO_'},
                  'ECO': {'prefix_prefix': 'ECO',
                          'prefix_reference': 'http://purl.obolibrary.org/obo/ECO_'},
                  'GO': {'prefix_prefix': 'GO',
                         'prefix_reference': 'http://purl.obolibrary.org/obo/GO_'},
                  'NCBITaxon': {'prefix_prefix': 'NCBITaxon',
                                'prefix_reference': 'http://purl.obolibrary.org/obo/NCBITaxon_'},
                  'OBAN': {'prefix_prefix': 'OBAN',
                           'prefix_reference': 'http://purl.org/oban/'},
                  'PMID': {'prefix_prefix': 'PMID',
                           'prefix_reference': 'http://identifiers.org/pmid/'},
                  'RO': {'prefix_prefix': 'RO',
                         'prefix_reference': 'http://purl.obolibrary.org/obo/RO_'},
                  'UniProtKB': {'prefix_prefix': 'UniProtKB',
                                'prefix_reference': 'http://identifiers.org/uniprot/'},
                  'biolink': {'prefix_prefix': 'biolink',
                              'prefix_reference': 'https://w3id.org/biolink/vocab/'},
                  'dce': {'prefix_prefix': 'dce',
                          'prefix_reference': 'http://purl.org/dc/elements/1.1/'},
                  'dcterms': {'prefix_prefix': 'dcterms',
                              'prefix_reference': 'http://purl.org/dc/terms/'},
                  'gocam': {'prefix_prefix': 'gocam',
                            'prefix_reference': 'https://w3id.org/gocam/'},
                  'gomodel': {'prefix_prefix': 'gomodel',
                              'prefix_reference': 'http://model.geneontology.org/'},
                  'goshapes': {'prefix_prefix': 'goshapes',
                               'prefix_reference': 'http://purl.obolibrary.org/obo/go/shapes/'},
                  'lego': {'prefix_prefix': 'lego',
                           'prefix_reference': 'http://geneontology.org/lego/'},
                  'linkml': {'prefix_prefix': 'linkml',
                             'prefix_reference': 'https://w3id.org/linkml/'},
                  'oio': {'prefix_prefix': 'oio',
                          'prefix_reference': 'http://www.geneontology.org/formats/oboInOwl#'},
                  'orcid': {'prefix_prefix': 'orcid',
                            'prefix_reference': 'https://orcid.org/'},
                  'pav': {'prefix_prefix': 'pav',
                          'prefix_reference': 'http://purl.org/pav/'}},
     'see_also': ['https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7012280/',
                  'https://docs.google.com/presentation/d/1ja0Vkw0AoENJ58emM77dGnqPtY1nfIJMeyVnObBxIxI/edit#slide=id.p8'],
     'source_file': 'src/gocam/schema/gocam.yaml'} )

class ModelStateEnum(str, Enum):
    """
    Status of a model
    """
    production = "production"
    development = "development"


class InformationBiomacromoleculeCategory(str, Enum):
    GeneOrReferenceProtein = "GeneOrReferenceProtein"
    ProteinIsoform = "ProteinIsoform"
    MacromolecularComplex = "MacromolecularComplex"
    Unknown = "Unknown"


class CausalPredicateEnum(str, Enum):
    causally_upstream_of_positive_effect = "causally upstream of, positive effect"
    causally_upstream_of_negative_effect = "causally upstream of, negative effect"
    causally_upstream_of = "causally upstream of"
    immediately_causally_upstream_of = "immediately causally upstream of"
    causally_upstream_of_or_within = "causally upstream of or within"
    causally_upstream_of_or_within_negative_effect = "causally upstream of or within, negative effect"
    causally_upstream_of_or_within_positive_effect = "causally upstream of or within, positive effect"
    regulates = "regulates"
    negatively_regulates = "negatively regulates"
    positively_regulates = "positively regulates"



class Model(ConfiguredBaseModel):
    """
    A model of a biological program consisting of a set of causally connected activities
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam'})

    id: str = Field(default=..., description="""The identifier of the model. Should be in gocam namespace.""", json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    title: Optional[str] = Field(default=None, description="""The human-readable descriptive title of the model""", json_schema_extra = { "linkml_meta": {'alias': 'title', 'domain_of': ['Model'], 'slot_uri': 'dct:title'} })
    taxon: Optional[str] = Field(default=None, description="""The primary taxon that the model is about""", json_schema_extra = { "linkml_meta": {'alias': 'taxon', 'domain_of': ['Model']} })
    status: Optional[ModelStateEnum] = Field(default=None, description="""The status of the model""", json_schema_extra = { "linkml_meta": {'alias': 'status',
         'aliases': ['model state'],
         'domain_of': ['Model'],
         'slot_uri': 'pav:status'} })
    comments: Optional[List[str]] = Field(default=None, description="""Comments about the model""", json_schema_extra = { "linkml_meta": {'alias': 'comments', 'domain_of': ['Model'], 'slot_uri': 'rdfs:comment'} })
    activities: Optional[List[Activity]] = Field(default=None, description="""All of the activities that are part of the model""", json_schema_extra = { "linkml_meta": {'alias': 'activities', 'domain_of': ['Model']} })
    objects: Optional[List[Union[Object,TermObject,PublicationObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""All of the objects that are part of the model. This includes terms as well as publications and database objects like gene. This is not strictly part of the data managed by the model, it is for convenience, and should be refreshed from outside.""", json_schema_extra = { "linkml_meta": {'alias': 'objects', 'domain_of': ['Model']} })
    provenances: Optional[List[ProvenanceInfo]] = Field(default=None, description="""Model-level provenance information""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class Activity(ConfiguredBaseModel):
    """
    An individual activity in a causal model, representing the individual molecular activity of a single gene product or complex in the context of a particular model
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'aliases': ['annoton'], 'from_schema': 'https://w3id.org/gocam'})

    id: str = Field(default=..., description="""Identifier of the activity unit. Should be in gocam namespace.""", json_schema_extra = { "linkml_meta": {'alias': 'id',
         'comments': ['Typically does not need to be exposed to end-user, this exists '
                      'to allow activity flows'],
         'domain_of': ['Model', 'Activity', 'Object']} })
    enabled_by: Optional[Union[EnabledByAssociation,EnabledByGeneProductAssociation,EnabledByProteinComplexAssociation]] = Field(default=None, description="""The gene product or complex that carries out the activity""", json_schema_extra = { "linkml_meta": {'alias': 'enabled_by', 'domain_of': ['Activity']} })
    molecular_function: Optional[MolecularFunctionAssociation] = Field(default=None, description="""The molecular function that is carried out by the gene product or complex""", json_schema_extra = { "linkml_meta": {'alias': 'molecular_function',
         'domain_of': ['Activity'],
         'todos': ['currently BP, CC etc are at the level of the activity, not the '
                   'MolecularFunctionAssociation']} })
    occurs_in: Optional[CellularAnatomicalEntityAssociation] = Field(default=None, description="""The cellular location in which the activity occurs""", json_schema_extra = { "linkml_meta": {'alias': 'occurs_in', 'domain_of': ['Activity']} })
    part_of: Optional[BiologicalProcessAssociation] = Field(default=None, description="""The larger biological process in which the activity is a part""", json_schema_extra = { "linkml_meta": {'alias': 'part_of',
         'domain_of': ['Activity',
                       'BiologicalProcessAssociation',
                       'CellularAnatomicalEntityAssociation',
                       'CellTypeAssociation',
                       'GrossAnatomyAssociation']} })
    has_input: Optional[List[MoleculeAssociation]] = Field(default=None, description="""The input molecules that are directly consumed by the activity""", json_schema_extra = { "linkml_meta": {'alias': 'has_input', 'domain_of': ['Activity']} })
    has_primary_input: Optional[MoleculeAssociation] = Field(default=None, description="""The primary input molecule that is directly consumed by the activity""", json_schema_extra = { "linkml_meta": {'alias': 'has_primary_input', 'domain_of': ['Activity']} })
    has_output: Optional[List[MoleculeAssociation]] = Field(default=None, description="""The output molecules that are directly produced by the activity""", json_schema_extra = { "linkml_meta": {'alias': 'has_output', 'domain_of': ['Activity']} })
    has_primary_output: Optional[MoleculeAssociation] = Field(default=None, description="""The primary output molecule that is directly produced by the activity""", json_schema_extra = { "linkml_meta": {'alias': 'has_primary_output', 'domain_of': ['Activity']} })
    causal_associations: Optional[List[CausalAssociation]] = Field(default=None, description="""The causal associations that connect this activity to other activities""", json_schema_extra = { "linkml_meta": {'alias': 'causal_associations', 'domain_of': ['Activity']} })
    provenances: Optional[List[ProvenanceInfo]] = Field(default=None, description="""Provenance information for the activity""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class EvidenceItem(ConfiguredBaseModel):
    """
    An individual piece of evidence that is associated with an assertion in a model
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam'})

    term: Optional[str] = Field(default=None, description="""The ECO term representing the type of evidence""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    reference: Optional[str] = Field(default=None, description="""The publication of reference that describes the evidence""", json_schema_extra = { "linkml_meta": {'alias': 'reference', 'domain_of': ['EvidenceItem']} })
    with_objects: Optional[List[str]] = Field(default=None, description="""Supporting database entities or terms""", json_schema_extra = { "linkml_meta": {'alias': 'with_objects',
         'aliases': ['with', 'with/from'],
         'domain_of': ['EvidenceItem']} })
    provenances: Optional[List[ProvenanceInfo]] = Field(default=None, description="""Provenance about the assertion, e.g. who made it""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class Association(ConfiguredBaseModel):
    """
    An abstract grouping for different kinds of evidence-associated provenance
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'abstract': True, 'from_schema': 'https://w3id.org/gocam'})

    type: Literal["Association"] = Field(default="Association", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[List[EvidenceItem]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[List[ProvenanceInfo]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class EnabledByAssociation(Association):
    """
    An association between an activity and the gene product or complex that carries it out
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'abstract': True, 'from_schema': 'https://w3id.org/gocam'})

    term: Optional[str] = Field(default=None, description="""The gene product or complex that carries out the activity""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    type: Literal["EnabledByAssociation"] = Field(default="EnabledByAssociation", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[List[EvidenceItem]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[List[ProvenanceInfo]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class EnabledByGeneProductAssociation(EnabledByAssociation):
    """
    An association between an activity and a gene product
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'slot_usage': {'term': {'name': 'term', 'range': 'GeneProductTermObject'}}})

    term: Optional[str] = Field(default=None, description="""The gene product or complex that carries out the activity""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    type: Literal["EnabledByGeneProductAssociation"] = Field(default="EnabledByGeneProductAssociation", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[List[EvidenceItem]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[List[ProvenanceInfo]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class EnabledByProteinComplexAssociation(EnabledByAssociation):
    """
    An association between an activity and a protein complex
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'slot_usage': {'term': {'name': 'term', 'range': 'ProteinComplexTermObject'}}})

    members: Optional[List[str]] = Field(default=None, description="""The gene products that are part of the complex""", json_schema_extra = { "linkml_meta": {'alias': 'members', 'domain_of': ['EnabledByProteinComplexAssociation']} })
    term: Optional[str] = Field(default=None, description="""The gene product or complex that carries out the activity""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    type: Literal["EnabledByProteinComplexAssociation"] = Field(default="EnabledByProteinComplexAssociation", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[List[EvidenceItem]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[List[ProvenanceInfo]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class CausalAssociation(Association):
    """
    A causal association between two activities
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam'})

    predicate: Optional[str] = Field(default=None, description="""The RO relation that represents the type of relationship""", json_schema_extra = { "linkml_meta": {'alias': 'predicate', 'domain_of': ['CausalAssociation']} })
    downstream_activity: Optional[str] = Field(default=None, description="""The activity unit that is downstream of this one""", json_schema_extra = { "linkml_meta": {'alias': 'downstream_activity',
         'aliases': ['object'],
         'domain_of': ['CausalAssociation']} })
    type: Literal["CausalAssociation"] = Field(default="CausalAssociation", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[List[EvidenceItem]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[List[ProvenanceInfo]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class TermAssociation(Association):
    """
    An association between an activity and a term, potentially with extensions
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'abstract': True, 'from_schema': 'https://w3id.org/gocam'})

    term: Optional[str] = Field(default=None, description="""The ontology term that describes the nature of the association""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    type: Literal["TermAssociation"] = Field(default="TermAssociation", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[List[EvidenceItem]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[List[ProvenanceInfo]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class MolecularFunctionAssociation(TermAssociation):
    """
    An association between an activity and a molecular function term
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'slot_usage': {'term': {'name': 'term',
                                 'range': 'MolecularFunctionTermObject'}},
         'todos': ['account for non-MF activity types in Reactome']})

    term: Optional[str] = Field(default=None, description="""The ontology term that describes the nature of the association""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    type: Literal["MolecularFunctionAssociation"] = Field(default="MolecularFunctionAssociation", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[List[EvidenceItem]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[List[ProvenanceInfo]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class BiologicalProcessAssociation(TermAssociation):
    """
    An association between an activity and a biological process term
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'slot_usage': {'term': {'name': 'term',
                                 'range': 'BiologicalProcessTermObject'}}})

    happens_during: Optional[str] = Field(default=None, description="""Optional extension describing where the BP takes place""", json_schema_extra = { "linkml_meta": {'alias': 'happens_during', 'domain_of': ['BiologicalProcessAssociation']} })
    part_of: Optional[str] = Field(default=None, description="""Optional extension allowing hierarchical nesting of BPs""", json_schema_extra = { "linkml_meta": {'alias': 'part_of',
         'domain_of': ['Activity',
                       'BiologicalProcessAssociation',
                       'CellularAnatomicalEntityAssociation',
                       'CellTypeAssociation',
                       'GrossAnatomyAssociation']} })
    term: Optional[str] = Field(default=None, description="""The ontology term that describes the nature of the association""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    type: Literal["BiologicalProcessAssociation"] = Field(default="BiologicalProcessAssociation", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[List[EvidenceItem]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[List[ProvenanceInfo]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class CellularAnatomicalEntityAssociation(TermAssociation):
    """
    An association between an activity and a cellular anatomical entity term
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'slot_usage': {'term': {'name': 'term',
                                 'range': 'CellularAnatomicalEntityTermObject'}}})

    part_of: Optional[CellTypeAssociation] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'part_of',
         'domain_of': ['Activity',
                       'BiologicalProcessAssociation',
                       'CellularAnatomicalEntityAssociation',
                       'CellTypeAssociation',
                       'GrossAnatomyAssociation']} })
    term: Optional[str] = Field(default=None, description="""The ontology term that describes the nature of the association""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    type: Literal["CellularAnatomicalEntityAssociation"] = Field(default="CellularAnatomicalEntityAssociation", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[List[EvidenceItem]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[List[ProvenanceInfo]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class CellTypeAssociation(TermAssociation):
    """
    An association between an activity and a cell type term
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'slot_usage': {'term': {'name': 'term', 'range': 'CellTypeTermObject'}}})

    part_of: Optional[GrossAnatomyAssociation] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'part_of',
         'domain_of': ['Activity',
                       'BiologicalProcessAssociation',
                       'CellularAnatomicalEntityAssociation',
                       'CellTypeAssociation',
                       'GrossAnatomyAssociation']} })
    term: Optional[str] = Field(default=None, description="""The ontology term that describes the nature of the association""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    type: Literal["CellTypeAssociation"] = Field(default="CellTypeAssociation", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[List[EvidenceItem]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[List[ProvenanceInfo]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class GrossAnatomyAssociation(TermAssociation):
    """
    An association between an activity and a gross anatomical structure term
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'slot_usage': {'term': {'name': 'term',
                                 'range': 'GrossAnatomicalStructureTermObject'}}})

    part_of: Optional[GrossAnatomyAssociation] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'part_of',
         'domain_of': ['Activity',
                       'BiologicalProcessAssociation',
                       'CellularAnatomicalEntityAssociation',
                       'CellTypeAssociation',
                       'GrossAnatomyAssociation']} })
    term: Optional[str] = Field(default=None, description="""The ontology term that describes the nature of the association""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    type: Literal["GrossAnatomyAssociation"] = Field(default="GrossAnatomyAssociation", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[List[EvidenceItem]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[List[ProvenanceInfo]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class MoleculeAssociation(TermAssociation):
    """
    An association between an activity and a molecule term
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'slot_usage': {'term': {'name': 'term', 'range': 'MoleculeTermObject'}}})

    term: Optional[str] = Field(default=None, description="""The ontology term that describes the nature of the association""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    type: Literal["MoleculeAssociation"] = Field(default="MoleculeAssociation", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[List[EvidenceItem]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[List[ProvenanceInfo]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class Object(ConfiguredBaseModel):
    """
    An abstract class for all identified objects in a model
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam'})

    id: str = Field(default=..., json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    label: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'label', 'domain_of': ['Object'], 'slot_uri': 'rdfs:label'} })
    type: Literal["https://w3id.org/gocam/Object","gocam:Object"] = Field(default="gocam:Object", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    obsolete: Optional[bool] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'obsolete', 'domain_of': ['Object']} })


class TermObject(Object):
    """
    An abstract class for all ontology term objects
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'abstract': True, 'from_schema': 'https://w3id.org/gocam'})

    id: str = Field(default=..., json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    label: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'label', 'domain_of': ['Object'], 'slot_uri': 'rdfs:label'} })
    type: Literal["https://w3id.org/gocam/TermObject","gocam:TermObject"] = Field(default="gocam:TermObject", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    obsolete: Optional[bool] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'obsolete', 'domain_of': ['Object']} })


class PublicationObject(Object):
    """
    An object that represents a publication or other kind of reference
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'id_prefixes': ['PMID', 'GOREF', 'DOI']})

    abstract_text: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'abstract_text', 'domain_of': ['PublicationObject']} })
    full_text: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'full_text', 'domain_of': ['PublicationObject']} })
    id: str = Field(default=..., json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    label: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'label', 'domain_of': ['Object'], 'slot_uri': 'rdfs:label'} })
    type: Literal["https://w3id.org/gocam/PublicationObject","gocam:PublicationObject"] = Field(default="gocam:PublicationObject", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    obsolete: Optional[bool] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'obsolete', 'domain_of': ['Object']} })


class EvidenceTermObject(TermObject):
    """
    A term object that represents an evidence term from ECO
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam', 'id_prefixes': ['ECO']})

    id: str = Field(default=..., json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    label: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'label', 'domain_of': ['Object'], 'slot_uri': 'rdfs:label'} })
    type: Literal["https://w3id.org/gocam/EvidenceTermObject","gocam:EvidenceTermObject"] = Field(default="gocam:EvidenceTermObject", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    obsolete: Optional[bool] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'obsolete', 'domain_of': ['Object']} })


class MolecularFunctionTermObject(TermObject):
    """
    A term object that represents a molecular function term from GO
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam', 'id_prefixes': ['GO']})

    id: str = Field(default=..., json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    label: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'label', 'domain_of': ['Object'], 'slot_uri': 'rdfs:label'} })
    type: Literal["https://w3id.org/gocam/MolecularFunctionTermObject","gocam:MolecularFunctionTermObject"] = Field(default="gocam:MolecularFunctionTermObject", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    obsolete: Optional[bool] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'obsolete', 'domain_of': ['Object']} })


class BiologicalProcessTermObject(TermObject):
    """
    A term object that represents a biological process term from GO
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam', 'id_prefixes': ['GO']})

    id: str = Field(default=..., json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    label: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'label', 'domain_of': ['Object'], 'slot_uri': 'rdfs:label'} })
    type: Literal["https://w3id.org/gocam/BiologicalProcessTermObject","gocam:BiologicalProcessTermObject"] = Field(default="gocam:BiologicalProcessTermObject", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    obsolete: Optional[bool] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'obsolete', 'domain_of': ['Object']} })


class CellularAnatomicalEntityTermObject(TermObject):
    """
    A term object that represents a cellular anatomical entity term from GO
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam', 'id_prefixes': ['GO']})

    id: str = Field(default=..., json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    label: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'label', 'domain_of': ['Object'], 'slot_uri': 'rdfs:label'} })
    type: Literal["https://w3id.org/gocam/CellularAnatomicalEntityTermObject","gocam:CellularAnatomicalEntityTermObject"] = Field(default="gocam:CellularAnatomicalEntityTermObject", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    obsolete: Optional[bool] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'obsolete', 'domain_of': ['Object']} })


class MoleculeTermObject(TermObject):
    """
    A term object that represents a molecule term from CHEBI or UniProtKB
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam', 'id_prefixes': ['CHEBI', 'UniProtKB']})

    id: str = Field(default=..., json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    label: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'label', 'domain_of': ['Object'], 'slot_uri': 'rdfs:label'} })
    type: Literal["https://w3id.org/gocam/MoleculeTermObject","gocam:MoleculeTermObject"] = Field(default="gocam:MoleculeTermObject", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    obsolete: Optional[bool] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'obsolete', 'domain_of': ['Object']} })


class CellTypeTermObject(TermObject):
    """
    A term object that represents a cell type term from CL
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'id_prefixes': ['CL', 'PO', 'FAO', 'DDANAT']})

    id: str = Field(default=..., json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    label: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'label', 'domain_of': ['Object'], 'slot_uri': 'rdfs:label'} })
    type: Literal["https://w3id.org/gocam/CellTypeTermObject","gocam:CellTypeTermObject"] = Field(default="gocam:CellTypeTermObject", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    obsolete: Optional[bool] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'obsolete', 'domain_of': ['Object']} })


class GrossAnatomicalStructureTermObject(TermObject):
    """
    A term object that represents a gross anatomical structure term from UBERON
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'id_prefixes': ['UBERON', 'PO', 'FAO', 'DDANAT']})

    id: str = Field(default=..., json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    label: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'label', 'domain_of': ['Object'], 'slot_uri': 'rdfs:label'} })
    type: Literal["https://w3id.org/gocam/GrossAnatomicalStructureTermObject","gocam:GrossAnatomicalStructureTermObject"] = Field(default="gocam:GrossAnatomicalStructureTermObject", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    obsolete: Optional[bool] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'obsolete', 'domain_of': ['Object']} })


class PhaseTermObject(TermObject):
    """
    A term object that represents a phase term from GO or UBERON
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam', 'id_prefixes': ['GO', 'UBERON', 'PO']})

    id: str = Field(default=..., json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    label: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'label', 'domain_of': ['Object'], 'slot_uri': 'rdfs:label'} })
    type: Literal["https://w3id.org/gocam/PhaseTermObject","gocam:PhaseTermObject"] = Field(default="gocam:PhaseTermObject", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    obsolete: Optional[bool] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'obsolete', 'domain_of': ['Object']} })


class InformationBiomacromoleculeTermObject(TermObject):
    """
    An abstract class for all information biomacromolecule term objects
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'abstract': True, 'from_schema': 'https://w3id.org/gocam'})

    id: str = Field(default=..., json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    label: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'label', 'domain_of': ['Object'], 'slot_uri': 'rdfs:label'} })
    type: Literal["https://w3id.org/gocam/InformationBiomacromoleculeTermObject","gocam:InformationBiomacromoleculeTermObject"] = Field(default="gocam:InformationBiomacromoleculeTermObject", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    obsolete: Optional[bool] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'obsolete', 'domain_of': ['Object']} })


class GeneProductTermObject(InformationBiomacromoleculeTermObject):
    """
    A term object that represents a gene product term from GO or UniProtKB
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam'})

    id: str = Field(default=..., json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    label: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'label', 'domain_of': ['Object'], 'slot_uri': 'rdfs:label'} })
    type: Literal["https://w3id.org/gocam/GeneProductTermObject","gocam:GeneProductTermObject"] = Field(default="gocam:GeneProductTermObject", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    obsolete: Optional[bool] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'obsolete', 'domain_of': ['Object']} })


class ProteinComplexTermObject(InformationBiomacromoleculeTermObject):
    """
    A term object that represents a protein complex term from GO
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam'})

    id: str = Field(default=..., json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    label: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'label', 'domain_of': ['Object'], 'slot_uri': 'rdfs:label'} })
    type: Literal["https://w3id.org/gocam/ProteinComplexTermObject","gocam:ProteinComplexTermObject"] = Field(default="gocam:ProteinComplexTermObject", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    obsolete: Optional[bool] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'obsolete', 'domain_of': ['Object']} })


class TaxonTermObject(TermObject):
    """
    A term object that represents a taxon term from NCBITaxon
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam', 'id_prefixes': ['NCBITaxon']})

    id: str = Field(default=..., json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    label: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'label', 'domain_of': ['Object'], 'slot_uri': 'rdfs:label'} })
    type: Literal["https://w3id.org/gocam/TaxonTermObject","gocam:TaxonTermObject"] = Field(default="gocam:TaxonTermObject", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    obsolete: Optional[bool] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'obsolete', 'domain_of': ['Object']} })


class PredicateTermObject(TermObject):
    """
    A term object that represents a taxon term from NCBITaxon
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam', 'id_prefixes': ['RO']})

    id: str = Field(default=..., json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    label: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'label', 'domain_of': ['Object'], 'slot_uri': 'rdfs:label'} })
    type: Literal["https://w3id.org/gocam/PredicateTermObject","gocam:PredicateTermObject"] = Field(default="gocam:PredicateTermObject", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    obsolete: Optional[bool] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'obsolete', 'domain_of': ['Object']} })


class ProvenanceInfo(ConfiguredBaseModel):
    """
    Provenance information for an object
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam'})

    contributor: Optional[List[str]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'contributor',
         'domain_of': ['ProvenanceInfo'],
         'slot_uri': 'dct:contributor'} })
    created: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'created', 'domain_of': ['ProvenanceInfo'], 'slot_uri': 'dct:created'} })
    date: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'date', 'domain_of': ['ProvenanceInfo'], 'slot_uri': 'dct:date'} })
    provided_by: Optional[List[str]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'provided_by',
         'domain_of': ['ProvenanceInfo'],
         'slot_uri': 'pav:providedBy'} })


# Model rebuild
# see https://pydantic-docs.helpmanual.io/usage/models/#rebuilding-a-model
Model.model_rebuild()
Activity.model_rebuild()
EvidenceItem.model_rebuild()
Association.model_rebuild()
EnabledByAssociation.model_rebuild()
EnabledByGeneProductAssociation.model_rebuild()
EnabledByProteinComplexAssociation.model_rebuild()
CausalAssociation.model_rebuild()
TermAssociation.model_rebuild()
MolecularFunctionAssociation.model_rebuild()
BiologicalProcessAssociation.model_rebuild()
CellularAnatomicalEntityAssociation.model_rebuild()
CellTypeAssociation.model_rebuild()
GrossAnatomyAssociation.model_rebuild()
MoleculeAssociation.model_rebuild()
Object.model_rebuild()
TermObject.model_rebuild()
PublicationObject.model_rebuild()
EvidenceTermObject.model_rebuild()
MolecularFunctionTermObject.model_rebuild()
BiologicalProcessTermObject.model_rebuild()
CellularAnatomicalEntityTermObject.model_rebuild()
MoleculeTermObject.model_rebuild()
CellTypeTermObject.model_rebuild()
GrossAnatomicalStructureTermObject.model_rebuild()
PhaseTermObject.model_rebuild()
InformationBiomacromoleculeTermObject.model_rebuild()
GeneProductTermObject.model_rebuild()
ProteinComplexTermObject.model_rebuild()
TaxonTermObject.model_rebuild()
PredicateTermObject.model_rebuild()
ProvenanceInfo.model_rebuild()

