from __future__ import annotations 
from datetime import (
    datetime,
    date
)
from decimal import Decimal 
from enum import Enum 
import re
import sys
from typing import (
    Any,
    List,
    Literal,
    Dict,
    Optional,
    Union
)
from pydantic.version import VERSION  as PYDANTIC_VERSION 
if int(PYDANTIC_VERSION[0])>=2:
    from pydantic import (
        BaseModel,
        ConfigDict,
        Field,
        field_validator
    )
else:
    from pydantic import (
        BaseModel,
        Field,
        validator
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
    id: str = Field(...)
    title: Optional[str] = Field(None)
    taxon: Optional[str] = Field(None)
    status: Optional[ModelStateEnum] = Field(None)
    comments: Optional[List[str]] = Field(default_factory=list)
    activities: Optional[List[Activity]] = Field(default_factory=list)
    objects: Optional[List[Union[Object,TermObject,PublicationObject,EvidenceTermObject,MolecularFunctionTermObject,BiologicalProcessTermObject,CellularAnatomicalEntityTermObject,MoleculeTermObject,CellTypeTermObject,GrossAnatomicalStructureTermObject,PhaseTermObject,InformationBiomacromoleculeTermObject,TaxonTermObject,PredicateTermObject,GeneProductTermObject,ProteinComplexTermObject]]] = Field(default_factory=list)
    provenances: Optional[List[ProvenanceInfo]] = Field(default_factory=list)


class Activity(ConfiguredBaseModel):
    """
    An individual activity in a causal model, representing the individual molecular activity of a single gene product or complex
    """
    id: str = Field(...)
    enabled_by: Optional[str] = Field(None, description="""The gene product or complex that carries out the activity""")
    molecular_function: Optional[MolecularFunctionAssociation] = Field(None, description="""The molecular function that is carried out by the gene product or complex""")
    occurs_in: Optional[CellularAnatomicalEntityAssociation] = Field(None, description="""The cellular location in which the activity occurs""")
    part_of: Optional[BiologicalProcessAssociation] = Field(None, description="""The larger biological process in which the activity is a part""")
    has_direct_input: Optional[MoleculeAssociation] = Field(None, description="""The input molecules that are directly consumed by the activity""")
    causal_associations: Optional[List[CausalAssociation]] = Field(default_factory=list, description="""The causal associations that connect this activity to other activities""")
    provenances: Optional[List[ProvenanceInfo]] = Field(default_factory=list, description="""Provenance information for the activity""")


class EvidenceItem(ConfiguredBaseModel):
    term: Optional[str] = Field(None)
    reference: Optional[str] = Field(None)
    with_objects: Optional[List[str]] = Field(default_factory=list)
    provenances: Optional[List[ProvenanceInfo]] = Field(default_factory=list)


class Association(ConfiguredBaseModel):
    """
    An abstract grouping for different kinds of evidence-associated provenance
    """
    type: Literal["Association"] = Field("Association")
    evidence: Optional[List[EvidenceItem]] = Field(default_factory=list)
    provenances: Optional[List[ProvenanceInfo]] = Field(default_factory=list)


class CausalAssociation(Association):
    """
    A causal association between two activities
    """
    predicate: Optional[str] = Field(None)
    downstream_activity: Optional[str] = Field(None)
    type: Literal["CausalAssociation"] = Field("CausalAssociation")
    evidence: Optional[List[EvidenceItem]] = Field(default_factory=list)
    provenances: Optional[List[ProvenanceInfo]] = Field(default_factory=list)


class TermAssociation(Association):
    """
    An association between an activity and a term, potentially with extensions
    """
    term: Optional[str] = Field(None)
    type: Literal["TermAssociation"] = Field("TermAssociation")
    evidence: Optional[List[EvidenceItem]] = Field(default_factory=list)
    provenances: Optional[List[ProvenanceInfo]] = Field(default_factory=list)


class MolecularFunctionAssociation(TermAssociation):
    """
    An association between an activity and a molecular function term
    """
    term: Optional[str] = Field(None)
    type: Literal["MolecularFunctionAssociation"] = Field("MolecularFunctionAssociation")
    evidence: Optional[List[EvidenceItem]] = Field(default_factory=list)
    provenances: Optional[List[ProvenanceInfo]] = Field(default_factory=list)


class BiologicalProcessAssociation(TermAssociation):
    """
    An association between an activity and a biological process term
    """
    term: Optional[str] = Field(None)
    type: Literal["BiologicalProcessAssociation"] = Field("BiologicalProcessAssociation")
    evidence: Optional[List[EvidenceItem]] = Field(default_factory=list)
    provenances: Optional[List[ProvenanceInfo]] = Field(default_factory=list)


class CellularAnatomicalEntityAssociation(TermAssociation):
    """
    An association between an activity and a cellular anatomical entity term
    """
    term: Optional[str] = Field(None)
    type: Literal["CellularAnatomicalEntityAssociation"] = Field("CellularAnatomicalEntityAssociation")
    evidence: Optional[List[EvidenceItem]] = Field(default_factory=list)
    provenances: Optional[List[ProvenanceInfo]] = Field(default_factory=list)


class MoleculeAssociation(TermAssociation):
    """
    An association between an activity and a molecule term
    """
    term: Optional[str] = Field(None)
    type: Literal["MoleculeAssociation"] = Field("MoleculeAssociation")
    evidence: Optional[List[EvidenceItem]] = Field(default_factory=list)
    provenances: Optional[List[ProvenanceInfo]] = Field(default_factory=list)


class Object(ConfiguredBaseModel):
    """
    An abstract class for all identified objects in a model
    """
    id: str = Field(...)
    label: Optional[str] = Field(None)
    type: Literal["https://w3id.org/gocam/Object","gocam:Object"] = Field("gocam:Object")
    obsolete: Optional[bool] = Field(None)


class TermObject(Object):
    """
    An abstract class for all ontology term objects
    """
    id: str = Field(...)
    label: Optional[str] = Field(None)
    type: Literal["https://w3id.org/gocam/TermObject","gocam:TermObject"] = Field("gocam:TermObject")
    obsolete: Optional[bool] = Field(None)


class PublicationObject(Object):
    """
    An object that represents a publication or other kind of reference
    """
    abstract_text: Optional[str] = Field(None)
    full_text: Optional[str] = Field(None)
    id: str = Field(...)
    label: Optional[str] = Field(None)
    type: Literal["https://w3id.org/gocam/PublicationObject","gocam:PublicationObject"] = Field("gocam:PublicationObject")
    obsolete: Optional[bool] = Field(None)


class EvidenceTermObject(TermObject):
    """
    A term object that represents an evidence term from ECO
    """
    id: str = Field(...)
    label: Optional[str] = Field(None)
    type: Literal["https://w3id.org/gocam/EvidenceTermObject","gocam:EvidenceTermObject"] = Field("gocam:EvidenceTermObject")
    obsolete: Optional[bool] = Field(None)


class MolecularFunctionTermObject(TermObject):
    """
    A term object that represents a molecular function term from GO
    """
    id: str = Field(...)
    label: Optional[str] = Field(None)
    type: Literal["https://w3id.org/gocam/MolecularFunctionTermObject","gocam:MolecularFunctionTermObject"] = Field("gocam:MolecularFunctionTermObject")
    obsolete: Optional[bool] = Field(None)


class BiologicalProcessTermObject(TermObject):
    """
    A termm object that represents a biological process term from GO
    """
    id: str = Field(...)
    label: Optional[str] = Field(None)
    type: Literal["https://w3id.org/gocam/BiologicalProcessTermObject","gocam:BiologicalProcessTermObject"] = Field("gocam:BiologicalProcessTermObject")
    obsolete: Optional[bool] = Field(None)


class CellularAnatomicalEntityTermObject(TermObject):
    """
    A term object that represents a cellular anatomical entity term from GO
    """
    id: str = Field(...)
    label: Optional[str] = Field(None)
    type: Literal["https://w3id.org/gocam/CellularAnatomicalEntityTermObject","gocam:CellularAnatomicalEntityTermObject"] = Field("gocam:CellularAnatomicalEntityTermObject")
    obsolete: Optional[bool] = Field(None)


class MoleculeTermObject(TermObject):
    """
    A term object that represents a molecule term from CHEBI or UniProtKB
    """
    id: str = Field(...)
    label: Optional[str] = Field(None)
    type: Literal["https://w3id.org/gocam/MoleculeTermObject","gocam:MoleculeTermObject"] = Field("gocam:MoleculeTermObject")
    obsolete: Optional[bool] = Field(None)


class CellTypeTermObject(TermObject):
    """
    A term object that represents a cell type term from CL
    """
    id: str = Field(...)
    label: Optional[str] = Field(None)
    type: Literal["https://w3id.org/gocam/CellTypeTermObject","gocam:CellTypeTermObject"] = Field("gocam:CellTypeTermObject")
    obsolete: Optional[bool] = Field(None)


class GrossAnatomicalStructureTermObject(TermObject):
    """
    A term object that represents a gross anatomical structure term from UBERON
    """
    id: str = Field(...)
    label: Optional[str] = Field(None)
    type: Literal["https://w3id.org/gocam/GrossAnatomicalStructureTermObject","gocam:GrossAnatomicalStructureTermObject"] = Field("gocam:GrossAnatomicalStructureTermObject")
    obsolete: Optional[bool] = Field(None)


class PhaseTermObject(TermObject):
    """
    A term object that represents a phase term from GO or UBERON
    """
    id: str = Field(...)
    label: Optional[str] = Field(None)
    type: Literal["https://w3id.org/gocam/PhaseTermObject","gocam:PhaseTermObject"] = Field("gocam:PhaseTermObject")
    obsolete: Optional[bool] = Field(None)


class InformationBiomacromoleculeTermObject(TermObject):
    """
    An abstract class for all information biomacromolecule term objects
    """
    id: str = Field(...)
    label: Optional[str] = Field(None)
    type: Literal["https://w3id.org/gocam/InformationBiomacromoleculeTermObject","gocam:InformationBiomacromoleculeTermObject"] = Field("gocam:InformationBiomacromoleculeTermObject")
    obsolete: Optional[bool] = Field(None)


class GeneProductTermObject(InformationBiomacromoleculeTermObject):
    """
    A term object that represents a gene product term from GO or UniProtKB
    """
    id: str = Field(...)
    label: Optional[str] = Field(None)
    type: Literal["https://w3id.org/gocam/GeneProductTermObject","gocam:GeneProductTermObject"] = Field("gocam:GeneProductTermObject")
    obsolete: Optional[bool] = Field(None)


class ProteinComplexTermObject(InformationBiomacromoleculeTermObject):
    """
    A term object that represents a protein complex term from GO
    """
    id: str = Field(...)
    label: Optional[str] = Field(None)
    type: Literal["https://w3id.org/gocam/ProteinComplexTermObject","gocam:ProteinComplexTermObject"] = Field("gocam:ProteinComplexTermObject")
    obsolete: Optional[bool] = Field(None)


class TaxonTermObject(TermObject):
    """
    A term object that represents a taxon term from NCBITaxon
    """
    id: str = Field(...)
    label: Optional[str] = Field(None)
    type: Literal["https://w3id.org/gocam/TaxonTermObject","gocam:TaxonTermObject"] = Field("gocam:TaxonTermObject")
    obsolete: Optional[bool] = Field(None)


class PredicateTermObject(TermObject):
    """
    A term object that represents a taxon term from NCBITaxon
    """
    id: str = Field(...)
    label: Optional[str] = Field(None)
    type: Literal["https://w3id.org/gocam/PredicateTermObject","gocam:PredicateTermObject"] = Field("gocam:PredicateTermObject")
    obsolete: Optional[bool] = Field(None)


class ProvenanceInfo(ConfiguredBaseModel):
    """
    Provenance information for an object
    """
    contributor: Optional[str] = Field(None)
    created: Optional[str] = Field(None)
    date: Optional[str] = Field(None)
    provided_by: Optional[str] = Field(None)


# Model rebuild
# see https://pydantic-docs.helpmanual.io/usage/models/#rebuilding-a-model
Model.model_rebuild()
Activity.model_rebuild()
EvidenceItem.model_rebuild()
Association.model_rebuild()
CausalAssociation.model_rebuild()
TermAssociation.model_rebuild()
MolecularFunctionAssociation.model_rebuild()
BiologicalProcessAssociation.model_rebuild()
CellularAnatomicalEntityAssociation.model_rebuild()
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
