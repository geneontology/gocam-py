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
    root: dict[str, Any] = {}
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
     'description': 'Gene Ontology Causal Activity Model (GO-CAM) Schema.\n'
                    '\n'
                    'This schema provides a way of representing causal pathway '
                    '[Models](Model.md). A model consists of a set of\n'
                    '[Activity](Activity.md) objects, where each activity object '
                    'represents the function of either an [individual\n'
                    'gene product](EnabledByGeneProductAssociation), a [protein '
                    'complex of gene products](EnabledByGeneProductAssociation),\n'
                    'or a set of possible gene products.\n'
                    '\n'
                    'Each [Models](Model.md) has associated metadata slots. Some '
                    'slots such as [id](id.md), [title](title.md),\n'
                    'and [status](status.md) are *required*.',
     'id': 'https://w3id.org/gocam',
     'imports': ['linkml:types'],
     'name': 'gocam',
     'prefixes': {'BFO': {'prefix_prefix': 'BFO',
                          'prefix_reference': 'http://purl.obolibrary.org/obo/BFO_'},
                  'CHEBI': {'prefix_prefix': 'CHEBI',
                            'prefix_reference': 'http://purl.obolibrary.org/obo/CHEBI_'},
                  'CL': {'prefix_prefix': 'CL',
                         'prefix_reference': 'http://purl.obolibrary.org/obo/CL_'},
                  'DDANAT': {'prefix_prefix': 'DDANAT',
                             'prefix_reference': 'http://purl.obolibrary.org/obo/DDANAT_'},
                  'DOI': {'prefix_prefix': 'DOI',
                          'prefix_reference': 'http://doi.org/'},
                  'ECO': {'prefix_prefix': 'ECO',
                          'prefix_reference': 'http://purl.obolibrary.org/obo/ECO_'},
                  'FAO': {'prefix_prefix': 'FAO',
                          'prefix_reference': 'http://purl.obolibrary.org/obo/FAO_'},
                  'GO': {'prefix_prefix': 'GO',
                         'prefix_reference': 'http://purl.obolibrary.org/obo/GO_'},
                  'GOREF': {'prefix_prefix': 'GOREF',
                            'prefix_reference': 'http://purl.obolibrary.org/obo/go/references/'},
                  'NCBITaxon': {'prefix_prefix': 'NCBITaxon',
                                'prefix_reference': 'http://purl.obolibrary.org/obo/NCBITaxon_'},
                  'OBAN': {'prefix_prefix': 'OBAN',
                           'prefix_reference': 'http://purl.org/oban/'},
                  'PMID': {'prefix_prefix': 'PMID',
                           'prefix_reference': 'http://identifiers.org/pubmed/'},
                  'PO': {'prefix_prefix': 'PO',
                         'prefix_reference': 'http://purl.obolibrary.org/obo/PO_'},
                  'RHEA': {'prefix_prefix': 'RHEA',
                           'prefix_reference': 'http://rdf.rhea-db.org/'},
                  'RO': {'prefix_prefix': 'RO',
                         'prefix_reference': 'http://purl.obolibrary.org/obo/RO_'},
                  'UBERON': {'prefix_prefix': 'UBERON',
                             'prefix_reference': 'http://purl.obolibrary.org/obo/UBERON_'},
                  'UniProtKB': {'prefix_prefix': 'UniProtKB',
                                'prefix_reference': 'http://purl.uniprot.org/uniprot/'},
                  'biolink': {'prefix_prefix': 'biolink',
                              'prefix_reference': 'https://w3id.org/biolink/vocab/'},
                  'dce': {'prefix_prefix': 'dce',
                          'prefix_reference': 'http://purl.org/dc/elements/1.1/'},
                  'dct': {'prefix_prefix': 'dct',
                          'prefix_reference': 'http://purl.org/dc/terms/'},
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
                          'prefix_reference': 'http://purl.org/pav/'},
                  'rdfs': {'prefix_prefix': 'rdfs',
                           'prefix_reference': 'http://www.w3.org/2000/01/rdf-schema#'}},
     'see_also': ['https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7012280/',
                  'https://docs.google.com/presentation/d/1ja0Vkw0AoENJ58emM77dGnqPtY1nfIJMeyVnObBxIxI/edit#slide=id.p8'],
     'source_file': 'src/gocam/schema/gocam.yaml'} )

class ModelStateEnum(str, Enum):
    """
    A term describing where the model is in the development life cycle.
    """
    development = "development"
    """
    Used when the curator is still working on the model. Edits are still being made, and the information in the model is not yet guaranteed to be accurate or complete. The model should not be displayed in end-user facing websites, unless it is made clear that the model is a work in progress.
    """
    production = "production"
    """
    Used when the curator has declared the model is ready for public consumption. Edits might still be performed on the model in future, but the information in the model is believed to be both accurate and reasonably complete. The model may be displayed in public websites.
    """
    delete = "delete"
    """
    When the curator has marked for future deletion.
    """
    review = "review"
    """
    The model has been marked for curator review.
    """
    internal_test = "internal_test"
    """
    The model is not intended for use public use; it is likely to be used for internal testing.
    """
    closed = "closed"
    """
    TBD
    """


class InformationBiomacromoleculeCategory(str, Enum):
    """
    A term describing the type of the enabler of an activity.
    """
    GeneOrReferenceProtein = "GeneOrReferenceProtein"
    ProteinIsoform = "ProteinIsoform"
    MacromolecularComplex = "MacromolecularComplex"
    Unknown = "Unknown"


class CausalPredicateEnum(str, Enum):
    """
    A term describing the causal relationship between two activities. All terms are drawn from the "causally upstream or within" (RO:0002418) branch of the Relation Ontology (RO).
    """
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
    provides_input_for = "provides input for"
    removes_input_for = "removes input for"


class EvidenceCodeEnum(str):
    """
    A term from the subset of ECO that maps up to a GAF evidence code
    """
    pass


class CellularAnatomicalEntityEnum(str):
    """
    A term from the subset of the cellular anatomical entity branch of GO CC
    """
    pass


class PhaseEnum(str):
    """
    A term from either the phase branch of GO or the phase branch of an anatomy ontology
    """
    pass



class Model(ConfiguredBaseModel):
    """
    A model of a biological program consisting of a set of causally connected activities.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'rules': [{'postconditions': {'slot_conditions': {'activities': {'name': 'activities',
                                                                          'required': True}}},
                    'preconditions': {'slot_conditions': {'status': {'equals_string': 'production',
                                                                     'name': 'status'}}},
                    'title': 'Production rules must have at least one activity'}]})

    id: str = Field(default=..., description="""The identifier of the model. Should be in gocam namespace.""", json_schema_extra = { "linkml_meta": {'alias': 'id', 'domain_of': ['Model', 'Activity', 'Object']} })
    title: str = Field(default=..., description="""The human-readable descriptive title of the model""", json_schema_extra = { "linkml_meta": {'alias': 'title', 'domain_of': ['Model'], 'slot_uri': 'dct:title'} })
    taxon: Optional[str] = Field(default=None, description="""The primary taxon that the model is about""", json_schema_extra = { "linkml_meta": {'alias': 'taxon', 'domain_of': ['Model']} })
    additional_taxa: Optional[list[str]] = Field(default=None, description="""Additional taxa that the model is about""", json_schema_extra = { "linkml_meta": {'alias': 'additional_taxa', 'domain_of': ['Model']} })
    status: Optional[ModelStateEnum] = Field(default=None, description="""The status of the model in terms of its progression along the developmental lifecycle""", json_schema_extra = { "linkml_meta": {'alias': 'status',
         'aliases': ['model state'],
         'domain_of': ['Model'],
         'slot_uri': 'pav:status'} })
    date_modified: Optional[str] = Field(default=None, description="""The date that the model was last modified""", json_schema_extra = { "linkml_meta": {'alias': 'date_modified', 'domain_of': ['Model'], 'slot_uri': 'dct:date'} })
    comments: Optional[list[str]] = Field(default=None, description="""Curator-provided comments about the model""", json_schema_extra = { "linkml_meta": {'alias': 'comments', 'domain_of': ['Model'], 'slot_uri': 'rdfs:comment'} })
    activities: Optional[list[Activity]] = Field(default=None, description="""All of the activities that are part of the model""", json_schema_extra = { "linkml_meta": {'alias': 'activities',
         'comments': ['this slot is conditionally required. It is optional for models '
                      'in development state (because a curator may need to instantiate '
                      'a Model before populating it with activities), but is required '
                      'for production models. See the associated rule.'],
         'domain_of': ['Model']} })
    objects: Optional[list[Union[Object,TermObject,PublicationObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""All of the objects that are part of the model. This includes terms as well as publications and database objects like gene. This is not strictly part of the data managed by the model, it is for convenience, and should be refreshed from outside.""", json_schema_extra = { "linkml_meta": {'alias': 'objects', 'domain_of': ['Model']} })
    provenances: Optional[list[ProvenanceInfo]] = Field(default=None, description="""Model-level provenance information""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })
    query_index: Optional[QueryIndex] = Field(default=None, description="""An optional object that contains the results of indexing a model with various summary statistics and retrieval indices.""", json_schema_extra = { "linkml_meta": {'alias': 'query_index',
         'comments': ['This is typically not populated in the primary transactional '
                      'store (OLTP processing), because the values will be redundant '
                      'with the primary edited components of the model. It is intended '
                      'to be populated in batch *after* editing, and then used for '
                      'generating reports, or for indexing in web applications.'],
         'domain_of': ['Model']} })


class Activity(ConfiguredBaseModel):
    """
    An individual activity in a causal model, representing the individual molecular activity of a single gene product or complex in the context of a particular model
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'aliases': ['annoton'], 'from_schema': 'https://w3id.org/gocam'})

    id: str = Field(default=..., description="""Identifier of the activity unit. Should be in gocam namespace.""", json_schema_extra = { "linkml_meta": {'alias': 'id',
         'comments': ['Typically does not need to be exposed to end-user, this exists '
                      'to allow activity flows'],
         'domain_of': ['Model', 'Activity', 'Object'],
         'examples': [{'description': 'A model representing tRNA repair and recycling',
                       'value': 'gomodel:63f809ec00000701'}],
         'id_prefixes': ['gocam']} })
    enabled_by: Optional[Union[EnabledByAssociation,EnabledByGeneProductAssociation,EnabledByProteinComplexAssociation]] = Field(default=None, description="""The gene product or complex that carries out the activity""", json_schema_extra = { "linkml_meta": {'alias': 'enabled_by', 'domain_of': ['Activity'], 'recommended': True} })
    molecular_function: Optional[MolecularFunctionAssociation] = Field(default=None, description="""The molecular function that is carried out by the gene product or complex""", json_schema_extra = { "linkml_meta": {'alias': 'molecular_function',
         'domain_of': ['Activity'],
         'recommended': True,
         'todos': ['currently BP, CC etc are at the level of the activity, not the '
                   'MolecularFunctionAssociation']} })
    occurs_in: Optional[CellularAnatomicalEntityAssociation] = Field(default=None, description="""The cellular location in which the activity occurs""", json_schema_extra = { "linkml_meta": {'alias': 'occurs_in', 'domain_of': ['Activity'], 'recommended': True} })
    part_of: Optional[BiologicalProcessAssociation] = Field(default=None, description="""The larger biological process in which the activity is a part""", json_schema_extra = { "linkml_meta": {'alias': 'part_of',
         'domain_of': ['Activity',
                       'BiologicalProcessAssociation',
                       'CellularAnatomicalEntityAssociation',
                       'CellTypeAssociation',
                       'GrossAnatomyAssociation'],
         'recommended': True} })
    has_input: Optional[list[MoleculeAssociation]] = Field(default=None, description="""The input molecules that are directly consumed by the activity""", json_schema_extra = { "linkml_meta": {'alias': 'has_input',
         'domain_of': ['Activity'],
         'todos': ['resolve has_input vs has_primary_input']} })
    has_primary_input: Optional[MoleculeAssociation] = Field(default=None, description="""The primary input molecule that is directly consumed by the activity""", json_schema_extra = { "linkml_meta": {'alias': 'has_primary_input',
         'domain_of': ['Activity'],
         'todos': ['resolve has_input vs has_primary_input']} })
    has_output: Optional[list[MoleculeAssociation]] = Field(default=None, description="""The output molecules that are directly produced by the activity""", json_schema_extra = { "linkml_meta": {'alias': 'has_output',
         'domain_of': ['Activity'],
         'todos': ['resolve has_output vs has_primary_output']} })
    has_primary_output: Optional[MoleculeAssociation] = Field(default=None, description="""The primary output molecule that is directly produced by the activity""", json_schema_extra = { "linkml_meta": {'alias': 'has_primary_output',
         'domain_of': ['Activity'],
         'todos': ['resolve has_output vs has_primary_output']} })
    causal_associations: Optional[list[CausalAssociation]] = Field(default=None, description="""The causal associations that flow out of this activity""", json_schema_extra = { "linkml_meta": {'alias': 'causal_associations',
         'comments': ['All activities in a model must be connected to at least one '
                      'other activity. If a an activity has no outgoing activities '
                      '(i.e the value of this slot is empty) then it is a terminal '
                      'activity in the model. If an activity has no incoming '
                      'activities, it is an initial activity.'],
         'domain_of': ['Activity']} })
    provenances: Optional[list[ProvenanceInfo]] = Field(default=None, description="""Provenance information for the activity""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class EvidenceItem(ConfiguredBaseModel):
    """
    An individual piece of evidence that is associated with an assertion in a model
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam'})

    term: Optional[str] = Field(default=None, description="""The ECO term representing the type of evidence""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'bindings': [{'binds_value_of': 'id',
                       'obligation_level': 'REQUIRED',
                       'range': 'EvidenceCodeEnum'}],
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation'],
         'examples': [{'description': 'direct assay evidence used in manual assertion '
                                      '(IDA)',
                       'value': 'ECO:0000314'}],
         'id_prefixes': ['ECO']} })
    reference: Optional[str] = Field(default=None, description="""The publication of reference that describes the evidence""", json_schema_extra = { "linkml_meta": {'alias': 'reference',
         'domain_of': ['EvidenceItem'],
         'examples': [{'value': 'PMID:32075755'}]} })
    with_objects: Optional[list[str]] = Field(default=None, description="""Supporting database entities or terms""", json_schema_extra = { "linkml_meta": {'alias': 'with_objects',
         'aliases': ['with', 'with/from'],
         'domain_of': ['EvidenceItem']} })
    provenances: Optional[list[ProvenanceInfo]] = Field(default=None, description="""Provenance about the assertion, e.g. who made it""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class Association(ConfiguredBaseModel):
    """
    An abstract grouping for different kinds of evidence-associated provenance
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'abstract': True, 'from_schema': 'https://w3id.org/gocam'})

    type: Literal["Association"] = Field(default="Association", description="""The type of association.""", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'comments': ['when instantiating Association objects in Python and other '
                      "languages, it isn't necessary to populate this, it is "
                      'auto-populated from the object class.'],
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[list[EvidenceItem]] = Field(default=None, description="""The set of evidence items that support the association.""", json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[list[ProvenanceInfo]] = Field(default=None, description="""The set of provenance objects that provide metadata about who made the association.""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class EnabledByAssociation(Association):
    """
    An association between an activity and the gene product or complex or set of potential gene products
      that carry out that activity.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'abstract': True,
         'comments': ['Note that this is an abstract class, and should ot be '
                      'instantiated directly, instead instantiate a subclass depending '
                      'on what kind of entity enables the association'],
         'from_schema': 'https://w3id.org/gocam'})

    term: Optional[str] = Field(default=None, description="""The gene product or complex that carries out the activity""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    type: Literal["EnabledByAssociation"] = Field(default="EnabledByAssociation", description="""The type of association.""", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'comments': ['when instantiating Association objects in Python and other '
                      "languages, it isn't necessary to populate this, it is "
                      'auto-populated from the object class.'],
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[list[EvidenceItem]] = Field(default=None, description="""The set of evidence items that support the association.""", json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[list[ProvenanceInfo]] = Field(default=None, description="""The set of provenance objects that provide metadata about who made the association.""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class EnabledByGeneProductAssociation(EnabledByAssociation):
    """
    An association between an activity and an individual gene product
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'slot_usage': {'term': {'comments': ['In the context of the GO workflow, the '
                                              'allowed values for this field come from '
                                              'the GPI file from an authoritative '
                                              'source. For example, the authoritative '
                                              'source for human is the EBI GOA group, '
                                              'and the GPI for this group consists of '
                                              'UniProtKB IDs (for proteins) and RNA '
                                              'Central IDs (for RNA gene products)',
                                              'A gene identifier may be provided as a '
                                              'value here (if the authoritative GPI '
                                              'allows it). Note that the '
                                              '*interpretation* of the gene ID in the '
                                              'context of a GO-CAM model is the '
                                              '(spliceform and proteoform agnostic) '
                                              '*product* of that gene.'],
                                 'description': 'A "term" that is an entity database '
                                                'object representing an individual '
                                                'gene product.',
                                 'examples': [{'description': 'The protein product of '
                                                              'the Homo sapiens TRNT1 '
                                                              'gene',
                                               'value': 'UniProtKB:Q96Q11'},
                                              {'description': 'An RNA product of this '
                                                              'RNA central gene',
                                               'value': 'RNAcentral:URS00026A1FBE_9606'}],
                                 'name': 'term',
                                 'range': 'GeneProductTermObject'}}})

    term: Optional[str] = Field(default=None, description="""A \"term\" that is an entity database object representing an individual gene product.""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'comments': ['In the context of the GO workflow, the allowed values for this '
                      'field come from the GPI file from an authoritative source. For '
                      'example, the authoritative source for human is the EBI GOA '
                      'group, and the GPI for this group consists of UniProtKB IDs '
                      '(for proteins) and RNA Central IDs (for RNA gene products)',
                      'A gene identifier may be provided as a value here (if the '
                      'authoritative GPI allows it). Note that the *interpretation* of '
                      'the gene ID in the context of a GO-CAM model is the (spliceform '
                      'and proteoform agnostic) *product* of that gene.'],
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation'],
         'examples': [{'description': 'The protein product of the Homo sapiens TRNT1 '
                                      'gene',
                       'value': 'UniProtKB:Q96Q11'},
                      {'description': 'An RNA product of this RNA central gene',
                       'value': 'RNAcentral:URS00026A1FBE_9606'}]} })
    type: Literal["EnabledByGeneProductAssociation"] = Field(default="EnabledByGeneProductAssociation", description="""The type of association.""", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'comments': ['when instantiating Association objects in Python and other '
                      "languages, it isn't necessary to populate this, it is "
                      'auto-populated from the object class.'],
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[list[EvidenceItem]] = Field(default=None, description="""The set of evidence items that support the association.""", json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[list[ProvenanceInfo]] = Field(default=None, description="""The set of provenance objects that provide metadata about who made the association.""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class EnabledByProteinComplexAssociation(EnabledByAssociation):
    """
    An association between an activity and a protein complex, where the complex carries out the activity. This should only be used when the activity cannot be attributed to an individual member of the complex, but instead the function is an emergent property of the complex.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'comments': ['Protein complexes can be specified either by *pre-composition* '
                      'or *post-composition*. For pre-composition, a species-specific '
                      'named protein complex (such as an entry in ComplexPortal) can '
                      'be specified, in which case the value of `members` is '
                      '*implicit*. For post-composition, the placeholder term '
                      '"GO:0032991" can be used, in which case `members` must be '
                      '*explicitly* specified. An intermediate case is when a named '
                      'class in GO that is a subclass of "GO:0032991" is used. In this '
                      'case, `members` should still be specified, as this may only be '
                      'partially specified by the GO class.'],
         'from_schema': 'https://w3id.org/gocam',
         'rules': [{'postconditions': {'slot_conditions': {'members': {'name': 'members',
                                                                       'required': True}}},
                    'preconditions': {'slot_conditions': {'term': {'equals_string': 'GO:0032991',
                                                                   'name': 'term'}}},
                    'title': 'members must be specified when the generic GO complex is '
                             'specified'}],
         'slot_usage': {'term': {'examples': [{'description': 'The generic GO entry '
                                                              'for a protein complex. '
                                                              'If this is the value of '
                                                              '`term`, then members '
                                                              '*must* be specified.',
                                               'value': 'GO:0032991'},
                                              {'description': 'The human Caspase-2 '
                                                              'complex',
                                               'value': 'ComplexPortal:CPX-969'}],
                                 'name': 'term',
                                 'range': 'ProteinComplexTermObject'}}})

    members: Optional[list[str]] = Field(default=None, description="""The gene products that are part of the complex""", json_schema_extra = { "linkml_meta": {'alias': 'members', 'domain_of': ['EnabledByProteinComplexAssociation']} })
    term: Optional[str] = Field(default=None, description="""The gene product or complex that carries out the activity""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation'],
         'examples': [{'description': 'The generic GO entry for a protein complex. If '
                                      'this is the value of `term`, then members '
                                      '*must* be specified.',
                       'value': 'GO:0032991'},
                      {'description': 'The human Caspase-2 complex',
                       'value': 'ComplexPortal:CPX-969'}]} })
    type: Literal["EnabledByProteinComplexAssociation"] = Field(default="EnabledByProteinComplexAssociation", description="""The type of association.""", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'comments': ['when instantiating Association objects in Python and other '
                      "languages, it isn't necessary to populate this, it is "
                      'auto-populated from the object class.'],
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[list[EvidenceItem]] = Field(default=None, description="""The set of evidence items that support the association.""", json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[list[ProvenanceInfo]] = Field(default=None, description="""The set of provenance objects that provide metadata about who made the association.""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
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
    type: Literal["CausalAssociation"] = Field(default="CausalAssociation", description="""The type of association.""", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'comments': ['when instantiating Association objects in Python and other '
                      "languages, it isn't necessary to populate this, it is "
                      'auto-populated from the object class.'],
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[list[EvidenceItem]] = Field(default=None, description="""The set of evidence items that support the association.""", json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[list[ProvenanceInfo]] = Field(default=None, description="""The set of provenance objects that provide metadata about who made the association.""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class TermAssociation(Association):
    """
    An association between an activity and a term, potentially with extensions. This is an abstract class for grouping purposes, it should not be directly instantiated, instead a subclass should be instantiated.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'abstract': True, 'from_schema': 'https://w3id.org/gocam'})

    term: Optional[str] = Field(default=None, description="""The ontology term that describes the nature of the association""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    type: Literal["TermAssociation"] = Field(default="TermAssociation", description="""The type of association.""", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'comments': ['when instantiating Association objects in Python and other '
                      "languages, it isn't necessary to populate this, it is "
                      'auto-populated from the object class.'],
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[list[EvidenceItem]] = Field(default=None, description="""The set of evidence items that support the association.""", json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[list[ProvenanceInfo]] = Field(default=None, description="""The set of provenance objects that provide metadata about who made the association.""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class MolecularFunctionAssociation(TermAssociation):
    """
    An association between an activity and a molecular function term
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'slot_usage': {'term': {'name': 'term',
                                 'range': 'MolecularFunctionTermObject'}},
         'todos': ['account for non-MF activity types in Reactome (MolecularEvent)']})

    term: Optional[str] = Field(default=None, description="""The ontology term that describes the nature of the association""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    type: Literal["MolecularFunctionAssociation"] = Field(default="MolecularFunctionAssociation", description="""The type of association.""", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'comments': ['when instantiating Association objects in Python and other '
                      "languages, it isn't necessary to populate this, it is "
                      'auto-populated from the object class.'],
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[list[EvidenceItem]] = Field(default=None, description="""The set of evidence items that support the association.""", json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[list[ProvenanceInfo]] = Field(default=None, description="""The set of provenance objects that provide metadata about who made the association.""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class BiologicalProcessAssociation(TermAssociation):
    """
    An association between an activity and a biological process term
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'slot_usage': {'term': {'name': 'term',
                                 'range': 'BiologicalProcessTermObject'}}})

    happens_during: Optional[str] = Field(default=None, description="""Optional extension describing where the BP takes place""", json_schema_extra = { "linkml_meta": {'alias': 'happens_during', 'domain_of': ['BiologicalProcessAssociation']} })
    part_of: Optional[BiologicalProcessAssociation] = Field(default=None, description="""Optional extension allowing hierarchical nesting of BPs""", json_schema_extra = { "linkml_meta": {'alias': 'part_of',
         'domain_of': ['Activity',
                       'BiologicalProcessAssociation',
                       'CellularAnatomicalEntityAssociation',
                       'CellTypeAssociation',
                       'GrossAnatomyAssociation']} })
    term: Optional[str] = Field(default=None, description="""The ontology term that describes the nature of the association""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    type: Literal["BiologicalProcessAssociation"] = Field(default="BiologicalProcessAssociation", description="""The type of association.""", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'comments': ['when instantiating Association objects in Python and other '
                      "languages, it isn't necessary to populate this, it is "
                      'auto-populated from the object class.'],
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[list[EvidenceItem]] = Field(default=None, description="""The set of evidence items that support the association.""", json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[list[ProvenanceInfo]] = Field(default=None, description="""The set of provenance objects that provide metadata about who made the association.""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class CellularAnatomicalEntityAssociation(TermAssociation):
    """
    An association between an activity and a cellular anatomical entity term
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'slot_usage': {'term': {'bindings': [{'binds_value_of': 'id',
                                               'obligation_level': 'REQUIRED',
                                               'range': 'CellularAnatomicalEntityEnum'}],
                                 'name': 'term',
                                 'range': 'CellularAnatomicalEntityTermObject'}}})

    part_of: Optional[CellTypeAssociation] = Field(default=None, description="""Optional extension allowing hierarchical nesting of CCs""", json_schema_extra = { "linkml_meta": {'alias': 'part_of',
         'domain_of': ['Activity',
                       'BiologicalProcessAssociation',
                       'CellularAnatomicalEntityAssociation',
                       'CellTypeAssociation',
                       'GrossAnatomyAssociation']} })
    term: Optional[str] = Field(default=None, description="""The ontology term that describes the nature of the association""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'bindings': [{'binds_value_of': 'id',
                       'obligation_level': 'REQUIRED',
                       'range': 'CellularAnatomicalEntityEnum'}],
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    type: Literal["CellularAnatomicalEntityAssociation"] = Field(default="CellularAnatomicalEntityAssociation", description="""The type of association.""", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'comments': ['when instantiating Association objects in Python and other '
                      "languages, it isn't necessary to populate this, it is "
                      'auto-populated from the object class.'],
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[list[EvidenceItem]] = Field(default=None, description="""The set of evidence items that support the association.""", json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[list[ProvenanceInfo]] = Field(default=None, description="""The set of provenance objects that provide metadata about who made the association.""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
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
    type: Literal["CellTypeAssociation"] = Field(default="CellTypeAssociation", description="""The type of association.""", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'comments': ['when instantiating Association objects in Python and other '
                      "languages, it isn't necessary to populate this, it is "
                      'auto-populated from the object class.'],
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[list[EvidenceItem]] = Field(default=None, description="""The set of evidence items that support the association.""", json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[list[ProvenanceInfo]] = Field(default=None, description="""The set of provenance objects that provide metadata about who made the association.""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
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
    type: Literal["GrossAnatomyAssociation"] = Field(default="GrossAnatomyAssociation", description="""The type of association.""", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'comments': ['when instantiating Association objects in Python and other '
                      "languages, it isn't necessary to populate this, it is "
                      'auto-populated from the object class.'],
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[list[EvidenceItem]] = Field(default=None, description="""The set of evidence items that support the association.""", json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[list[ProvenanceInfo]] = Field(default=None, description="""The set of provenance objects that provide metadata about who made the association.""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
         'domain_of': ['Model', 'Activity', 'EvidenceItem', 'Association']} })


class MoleculeAssociation(TermAssociation):
    """
    An association between an activity and a molecule term
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam',
         'slot_usage': {'term': {'name': 'term', 'range': 'MoleculeTermObject'}}})

    term: Optional[str] = Field(default=None, description="""The ontology term that describes the nature of the association""", json_schema_extra = { "linkml_meta": {'alias': 'term',
         'domain_of': ['EvidenceItem', 'EnabledByAssociation', 'TermAssociation']} })
    type: Literal["MoleculeAssociation"] = Field(default="MoleculeAssociation", description="""The type of association.""", json_schema_extra = { "linkml_meta": {'alias': 'type',
         'comments': ['when instantiating Association objects in Python and other '
                      "languages, it isn't necessary to populate this, it is "
                      'auto-populated from the object class.'],
         'designates_type': True,
         'domain_of': ['Association', 'Object']} })
    evidence: Optional[list[EvidenceItem]] = Field(default=None, description="""The set of evidence items that support the association.""", json_schema_extra = { "linkml_meta": {'alias': 'evidence', 'domain_of': ['Association']} })
    provenances: Optional[list[ProvenanceInfo]] = Field(default=None, description="""The set of provenance objects that provide metadata about who made the association.""", json_schema_extra = { "linkml_meta": {'alias': 'provenances',
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
    A term object that represents an evidence term from ECO. Only ECO terms that map up to a GO GAF evidence code should be used.
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

    contributor: Optional[list[str]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'contributor',
         'domain_of': ['ProvenanceInfo'],
         'slot_uri': 'dct:contributor'} })
    created: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'created', 'domain_of': ['ProvenanceInfo'], 'slot_uri': 'dct:created'} })
    date: Optional[str] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'date',
         'domain_of': ['ProvenanceInfo'],
         'slot_uri': 'dct:date',
         'todos': ['consider modeling as date rather than string']} })
    provided_by: Optional[list[str]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'provided_by',
         'domain_of': ['ProvenanceInfo'],
         'slot_uri': 'pav:providedBy'} })


class QueryIndex(ConfiguredBaseModel):
    """
    An index that is optionally placed on a model in order to support common query or index operations. Note that this index is not typically populated in the working transactional store for a model, it is derived via computation from core primary model information.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/gocam'})

    taxon_label: Optional[str] = Field(default=None, description="""The label of the primary taxon for the model""", json_schema_extra = { "linkml_meta": {'alias': 'taxon_label', 'domain_of': ['QueryIndex']} })
    number_of_activities: Optional[int] = Field(default=None, description="""The number of activities in a model.""", json_schema_extra = { "linkml_meta": {'alias': 'number_of_activities',
         'comments': ['this includes all activities, even those without an enabler.'],
         'domain_of': ['QueryIndex']} })
    number_of_enabled_by_terms: Optional[int] = Field(default=None, description="""The number of molecular entities or sets of entities in a model.""", json_schema_extra = { "linkml_meta": {'alias': 'number_of_enabled_by_terms', 'domain_of': ['QueryIndex']} })
    number_of_causal_associations: Optional[int] = Field(default=None, description="""Total number of causal association edges connecting activities in a model.""", json_schema_extra = { "linkml_meta": {'alias': 'number_of_causal_associations',
         'domain_of': ['QueryIndex'],
         'todos': ['decide what to do about "implicit" causal associations, i.e '
                   'provides_input_for']} })
    length_of_longest_causal_association_path: Optional[int] = Field(default=None, description="""The maximum number of hops along activities along the direction of causal flow in a model.""", json_schema_extra = { "linkml_meta": {'alias': 'length_of_longest_causal_association_path',
         'domain_of': ['QueryIndex']} })
    number_of_strongly_connected_components: Optional[int] = Field(default=None, description="""The number of distinct components that consist of activities that are connected (directly or indirectly) via causal connections. Most models will consist of a single SCC. Some models may consist of two or more \"islands\" where there is no connection from one island to another.""", json_schema_extra = { "linkml_meta": {'alias': 'number_of_strongly_connected_components',
         'domain_of': ['QueryIndex']} })
    flattened_references: Optional[list[PublicationObject]] = Field(default=None, description="""All publication objects from the model across different levels combined in one place""", json_schema_extra = { "linkml_meta": {'alias': 'flattened_references', 'domain_of': ['QueryIndex']} })
    model_activity_molecular_function_terms: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""All MF terms for all activities""", json_schema_extra = { "linkml_meta": {'alias': 'model_activity_molecular_function_terms',
         'domain_of': ['QueryIndex']} })
    model_activity_molecular_function_closure: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""The reflexive transitive closure of `model_activity_molecular_function_terms`, over the is_a relationship""", json_schema_extra = { "linkml_meta": {'alias': 'model_activity_molecular_function_closure',
         'annotations': {'closure_computed_over': {'tag': 'closure_computed_over',
                                                   'value': '[rdfs:subClassOf]'}},
         'domain_of': ['QueryIndex']} })
    model_activity_molecular_function_rollup: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""The rollup of `model_activity_molecular_function_closure` to a GO subset or slim.""", json_schema_extra = { "linkml_meta": {'alias': 'model_activity_molecular_function_rollup',
         'comments': ['added for completion but may not be useful in practice'],
         'domain_of': ['QueryIndex']} })
    model_activity_occurs_in_terms: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""All direct cellular component localization terms for all activities""", json_schema_extra = { "linkml_meta": {'alias': 'model_activity_occurs_in_terms', 'domain_of': ['QueryIndex']} })
    model_activity_occurs_in_closure: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""The reflexive transitive closure of `model_activity_occurs_in_terms`, over the is_a and part_of relationship type""", json_schema_extra = { "linkml_meta": {'alias': 'model_activity_occurs_in_closure',
         'annotations': {'closure_computed_over': {'tag': 'closure_computed_over',
                                                   'value': '[rdfs:subClassOf, '
                                                            'BFO:0000050]'}},
         'domain_of': ['QueryIndex']} })
    model_activity_occurs_in_rollup: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""The rollup of `model_activity_occurs_in_closure` to a GO subset or slim.""", json_schema_extra = { "linkml_meta": {'alias': 'model_activity_occurs_in_rollup', 'domain_of': ['QueryIndex']} })
    model_activity_enabled_by_terms: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""All direct enabler terms for all activities""", json_schema_extra = { "linkml_meta": {'alias': 'model_activity_enabled_by_terms', 'domain_of': ['QueryIndex']} })
    model_activity_enabled_by_closure: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""The reflexive transitive closure of `model_activity_enabled_by_terms`, over the is_a and part_of relationship types""", json_schema_extra = { "linkml_meta": {'alias': 'model_activity_enabled_by_closure',
         'annotations': {'closure_computed_over': {'tag': 'closure_computed_over',
                                                   'value': '[rdfs:subClassOf, '
                                                            'BFO:0000050]'}},
         'domain_of': ['QueryIndex']} })
    model_activity_enabled_by_rollup: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""The rollup of `model_activity_enabled_by_closure` to a GO subset or slim.""", json_schema_extra = { "linkml_meta": {'alias': 'model_activity_enabled_by_rollup',
         'comments': ['added for completion but may not be useful in practice'],
         'domain_of': ['QueryIndex']} })
    model_activity_enabled_by_genes: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""All direct enabler genes and genes that are part of enabler complexes for all activities""", json_schema_extra = { "linkml_meta": {'alias': 'model_activity_enabled_by_genes', 'domain_of': ['QueryIndex']} })
    model_activity_part_of_terms: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""All direct biological process terms for all activities""", json_schema_extra = { "linkml_meta": {'alias': 'model_activity_part_of_terms', 'domain_of': ['QueryIndex']} })
    model_activity_part_of_closure: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""The reflexive transitive closure of `model_activity_part_of_terms`, over the is_a and part_of relationship type""", json_schema_extra = { "linkml_meta": {'alias': 'model_activity_part_of_closure',
         'annotations': {'closure_computed_over': {'tag': 'closure_computed_over',
                                                   'value': '[rdfs:subClassOf, '
                                                            'BFO:0000050]'}},
         'domain_of': ['QueryIndex']} })
    model_activity_part_of_rollup: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""The rollup of `model_activity_part_of_closure` to a GO subset or slim.""", json_schema_extra = { "linkml_meta": {'alias': 'model_activity_part_of_rollup', 'domain_of': ['QueryIndex']} })
    model_activity_has_input_terms: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""All direct input terms for all activities""", json_schema_extra = { "linkml_meta": {'alias': 'model_activity_has_input_terms', 'domain_of': ['QueryIndex']} })
    model_activity_has_input_closure: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""The reflexive transitive closure of `model_activity_has_input_terms`, over the is_a relationship type""", json_schema_extra = { "linkml_meta": {'alias': 'model_activity_has_input_closure',
         'annotations': {'closure_computed_over': {'tag': 'closure_computed_over',
                                                   'value': '[rdfs:subClassOf]'}},
         'domain_of': ['QueryIndex']} })
    model_activity_has_input_rollup: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""The rollup of `model_activity_has_input_closure` to a GO subset or slim.""", json_schema_extra = { "linkml_meta": {'alias': 'model_activity_has_input_rollup', 'domain_of': ['QueryIndex']} })
    model_taxon: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""The primary taxon term for the model, over the NCBITaxon:subClassOf relationship type. This is used to determine the primary taxon that the model is relevant to.""", json_schema_extra = { "linkml_meta": {'alias': 'model_taxon', 'domain_of': ['QueryIndex']} })
    model_taxon_closure: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""The reflexive transitive closure of the taxon term for the model, over the NCBITaxon:subClassOf relationship type. This is used to determine the set of taxa that are relevant to the model.""", json_schema_extra = { "linkml_meta": {'alias': 'model_taxon_closure',
         'annotations': {'closure_computed_over': {'tag': 'closure_computed_over',
                                                   'value': '[rdfs:subClassOf]'}},
         'domain_of': ['QueryIndex']} })
    model_taxon_rollup: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, description="""The rollup of the taxon closure to a NCBITaxon subset or slim.""", json_schema_extra = { "linkml_meta": {'alias': 'model_taxon_rollup', 'domain_of': ['QueryIndex']} })
    annoton_terms: Optional[list[Union[TermObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default=None, json_schema_extra = { "linkml_meta": {'alias': 'annoton_terms', 'domain_of': ['QueryIndex']} })
    start_activities: Optional[list[str]] = Field(default=None, description="""The set of activities that are the starting points of the model, i.e. those that have no incoming causal associations.""", json_schema_extra = { "linkml_meta": {'alias': 'start_activities', 'domain_of': ['QueryIndex']} })
    end_activities: Optional[list[str]] = Field(default=None, description="""The set of activities that are the end points of the model, i.e. those that have no outgoing causal associations.""", json_schema_extra = { "linkml_meta": {'alias': 'end_activities', 'domain_of': ['QueryIndex']} })
    intermediate_activities: Optional[list[str]] = Field(default=None, description="""The set of activities that are neither start nor end activities, i.e. those that have both incoming and outgoing causal associations.""", json_schema_extra = { "linkml_meta": {'alias': 'intermediate_activities', 'domain_of': ['QueryIndex']} })
    singleton_activities: Optional[list[str]] = Field(default=None, description="""The set of activities that have no causal associations, i.e. those that are not connected to any other activity in the model.""", json_schema_extra = { "linkml_meta": {'alias': 'singleton_activities', 'domain_of': ['QueryIndex']} })
    number_of_start_activities: Optional[int] = Field(default=None, description="""The number of start activities in a model""", json_schema_extra = { "linkml_meta": {'alias': 'number_of_start_activities', 'domain_of': ['QueryIndex']} })
    number_of_end_activities: Optional[int] = Field(default=None, description="""The number of end activities in a model""", json_schema_extra = { "linkml_meta": {'alias': 'number_of_end_activities', 'domain_of': ['QueryIndex']} })
    number_of_intermediate_activities: Optional[int] = Field(default=None, description="""The number of intermediate activities in a model""", json_schema_extra = { "linkml_meta": {'alias': 'number_of_intermediate_activities', 'domain_of': ['QueryIndex']} })
    number_of_singleton_activities: Optional[int] = Field(default=None, description="""The number of singleton activities in a model""", json_schema_extra = { "linkml_meta": {'alias': 'number_of_singleton_activities', 'domain_of': ['QueryIndex']} })


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
QueryIndex.model_rebuild()

