-- # Class: "Model" Description: "A model of a biological program consisting of a set of causally connected activities"
--     * Slot: id Description: The identifier of the model. Should be in gocam namespace.
--     * Slot: title Description: The human-readable descriptive title of the model
--     * Slot: taxon Description: The primary taxon that the model is about
--     * Slot: status Description: The status of the model
-- # Class: "Activity" Description: "An individual activity in a causal model, representing the individual molecular activity of a single gene product or complex in the context of a particular model"
--     * Slot: id Description: Identifier of the activity unit. Should be in gocam namespace.
--     * Slot: Model_id Description: Autocreated FK slot
--     * Slot: enabled_by_id Description: The gene product or complex that carries out the activity
--     * Slot: molecular_function_id Description: The molecular function that is carried out by the gene product or complex
--     * Slot: occurs_in_id Description: The cellular location in which the activity occurs
--     * Slot: part_of_id Description: The larger biological process in which the activity is a part
--     * Slot: has_primary_input_id Description: The primary input molecule that is directly consumed by the activity
--     * Slot: has_primary_output_id Description: The primary output molecule that is directly produced by the activity
-- # Class: "EvidenceItem" Description: "An individual piece of evidence that is associated with an assertion in a model"
--     * Slot: id Description: 
--     * Slot: term Description: The ECO term representing the type of evidence
--     * Slot: reference Description: The publication of reference that describes the evidence
--     * Slot: Association_id Description: Autocreated FK slot
--     * Slot: EnabledByAssociation_id Description: Autocreated FK slot
--     * Slot: EnabledByGeneProductAssociation_id Description: Autocreated FK slot
--     * Slot: EnabledByProteinComplexAssociation_id Description: Autocreated FK slot
--     * Slot: CausalAssociation_id Description: Autocreated FK slot
--     * Slot: TermAssociation_id Description: Autocreated FK slot
--     * Slot: MolecularFunctionAssociation_id Description: Autocreated FK slot
--     * Slot: BiologicalProcessAssociation_id Description: Autocreated FK slot
--     * Slot: CellularAnatomicalEntityAssociation_id Description: Autocreated FK slot
--     * Slot: CellTypeAssociation_id Description: Autocreated FK slot
--     * Slot: GrossAnatomyAssociation_id Description: Autocreated FK slot
--     * Slot: MoleculeAssociation_id Description: Autocreated FK slot
-- # Class: "Association" Description: "An abstract grouping for different kinds of evidence-associated provenance"
--     * Slot: id Description: 
--     * Slot: type Description: 
-- # Class: "EnabledByAssociation" Description: "An association between an activity and the gene product or complex that carries it out"
--     * Slot: id Description: 
--     * Slot: term Description: The gene product or complex that carries out the activity
--     * Slot: type Description: 
-- # Class: "EnabledByGeneProductAssociation" Description: "An association between an activity and a gene product"
--     * Slot: id Description: 
--     * Slot: term Description: The gene product or complex that carries out the activity
--     * Slot: type Description: 
-- # Class: "EnabledByProteinComplexAssociation" Description: "An association between an activity and a protein complex"
--     * Slot: id Description: 
--     * Slot: term Description: The gene product or complex that carries out the activity
--     * Slot: type Description: 
-- # Class: "CausalAssociation" Description: "A causal association between two activities"
--     * Slot: id Description: 
--     * Slot: predicate Description: The RO relation that represents the type of relationship
--     * Slot: downstream_activity Description: The activity unit that is downstream of this one
--     * Slot: type Description: 
--     * Slot: Activity_id Description: Autocreated FK slot
-- # Class: "TermAssociation" Description: "An association between an activity and a term, potentially with extensions"
--     * Slot: id Description: 
--     * Slot: term Description: The ontology term that describes the nature of the association
--     * Slot: type Description: 
-- # Class: "MolecularFunctionAssociation" Description: "An association between an activity and a molecular function term"
--     * Slot: id Description: 
--     * Slot: term Description: The ontology term that describes the nature of the association
--     * Slot: type Description: 
-- # Class: "BiologicalProcessAssociation" Description: "An association between an activity and a biological process term"
--     * Slot: id Description: 
--     * Slot: happens_during Description: Optional extension describing where the BP takes place
--     * Slot: part_of Description: Optional extension allowing hierarchical nesting of BPs
--     * Slot: term Description: The ontology term that describes the nature of the association
--     * Slot: type Description: 
-- # Class: "CellularAnatomicalEntityAssociation" Description: "An association between an activity and a cellular anatomical entity term"
--     * Slot: id Description: 
--     * Slot: term Description: The ontology term that describes the nature of the association
--     * Slot: type Description: 
--     * Slot: part_of_id Description: 
-- # Class: "CellTypeAssociation" Description: "An association between an activity and a cell type term"
--     * Slot: id Description: 
--     * Slot: term Description: The ontology term that describes the nature of the association
--     * Slot: type Description: 
--     * Slot: part_of_id Description: 
-- # Class: "GrossAnatomyAssociation" Description: "An association between an activity and a gross anatomical structure term"
--     * Slot: id Description: 
--     * Slot: term Description: The ontology term that describes the nature of the association
--     * Slot: type Description: 
--     * Slot: part_of_id Description: 
-- # Class: "MoleculeAssociation" Description: "An association between an activity and a molecule term"
--     * Slot: id Description: 
--     * Slot: term Description: The ontology term that describes the nature of the association
--     * Slot: type Description: 
--     * Slot: Activity_id Description: Autocreated FK slot
-- # Class: "Object" Description: "An abstract class for all identified objects in a model"
--     * Slot: id Description: 
--     * Slot: label Description: 
--     * Slot: type Description: 
--     * Slot: obsolete Description: 
--     * Slot: Model_id Description: Autocreated FK slot
-- # Class: "TermObject" Description: "An abstract class for all ontology term objects"
--     * Slot: id Description: 
--     * Slot: label Description: 
--     * Slot: type Description: 
--     * Slot: obsolete Description: 
-- # Class: "PublicationObject" Description: "An object that represents a publication or other kind of reference"
--     * Slot: abstract_text Description: 
--     * Slot: full_text Description: 
--     * Slot: id Description: 
--     * Slot: label Description: 
--     * Slot: type Description: 
--     * Slot: obsolete Description: 
-- # Class: "EvidenceTermObject" Description: "A term object that represents an evidence term from ECO"
--     * Slot: id Description: 
--     * Slot: label Description: 
--     * Slot: type Description: 
--     * Slot: obsolete Description: 
-- # Class: "MolecularFunctionTermObject" Description: "A term object that represents a molecular function term from GO"
--     * Slot: id Description: 
--     * Slot: label Description: 
--     * Slot: type Description: 
--     * Slot: obsolete Description: 
-- # Class: "BiologicalProcessTermObject" Description: "A term object that represents a biological process term from GO"
--     * Slot: id Description: 
--     * Slot: label Description: 
--     * Slot: type Description: 
--     * Slot: obsolete Description: 
-- # Class: "CellularAnatomicalEntityTermObject" Description: "A term object that represents a cellular anatomical entity term from GO"
--     * Slot: id Description: 
--     * Slot: label Description: 
--     * Slot: type Description: 
--     * Slot: obsolete Description: 
-- # Class: "MoleculeTermObject" Description: "A term object that represents a molecule term from CHEBI or UniProtKB"
--     * Slot: id Description: 
--     * Slot: label Description: 
--     * Slot: type Description: 
--     * Slot: obsolete Description: 
-- # Class: "CellTypeTermObject" Description: "A term object that represents a cell type term from CL"
--     * Slot: id Description: 
--     * Slot: label Description: 
--     * Slot: type Description: 
--     * Slot: obsolete Description: 
-- # Class: "GrossAnatomicalStructureTermObject" Description: "A term object that represents a gross anatomical structure term from UBERON"
--     * Slot: id Description: 
--     * Slot: label Description: 
--     * Slot: type Description: 
--     * Slot: obsolete Description: 
-- # Class: "PhaseTermObject" Description: "A term object that represents a phase term from GO or UBERON"
--     * Slot: id Description: 
--     * Slot: label Description: 
--     * Slot: type Description: 
--     * Slot: obsolete Description: 
-- # Class: "InformationBiomacromoleculeTermObject" Description: "An abstract class for all information biomacromolecule term objects"
--     * Slot: id Description: 
--     * Slot: label Description: 
--     * Slot: type Description: 
--     * Slot: obsolete Description: 
-- # Class: "GeneProductTermObject" Description: "A term object that represents a gene product term from GO or UniProtKB"
--     * Slot: id Description: 
--     * Slot: label Description: 
--     * Slot: type Description: 
--     * Slot: obsolete Description: 
-- # Class: "ProteinComplexTermObject" Description: "A term object that represents a protein complex term from GO"
--     * Slot: id Description: 
--     * Slot: label Description: 
--     * Slot: type Description: 
--     * Slot: obsolete Description: 
-- # Class: "TaxonTermObject" Description: "A term object that represents a taxon term from NCBITaxon"
--     * Slot: id Description: 
--     * Slot: label Description: 
--     * Slot: type Description: 
--     * Slot: obsolete Description: 
-- # Class: "PredicateTermObject" Description: "A term object that represents a taxon term from NCBITaxon"
--     * Slot: id Description: 
--     * Slot: label Description: 
--     * Slot: type Description: 
--     * Slot: obsolete Description: 
-- # Class: "ProvenanceInfo" Description: "Provenance information for an object"
--     * Slot: id Description: 
--     * Slot: contributor Description: 
--     * Slot: created Description: 
--     * Slot: date Description: 
--     * Slot: provided_by Description: 
--     * Slot: Model_id Description: Autocreated FK slot
--     * Slot: Activity_id Description: Autocreated FK slot
--     * Slot: EvidenceItem_id Description: Autocreated FK slot
--     * Slot: Association_id Description: Autocreated FK slot
--     * Slot: EnabledByAssociation_id Description: Autocreated FK slot
--     * Slot: EnabledByGeneProductAssociation_id Description: Autocreated FK slot
--     * Slot: EnabledByProteinComplexAssociation_id Description: Autocreated FK slot
--     * Slot: CausalAssociation_id Description: Autocreated FK slot
--     * Slot: TermAssociation_id Description: Autocreated FK slot
--     * Slot: MolecularFunctionAssociation_id Description: Autocreated FK slot
--     * Slot: BiologicalProcessAssociation_id Description: Autocreated FK slot
--     * Slot: CellularAnatomicalEntityAssociation_id Description: Autocreated FK slot
--     * Slot: CellTypeAssociation_id Description: Autocreated FK slot
--     * Slot: GrossAnatomyAssociation_id Description: Autocreated FK slot
--     * Slot: MoleculeAssociation_id Description: Autocreated FK slot
-- # Class: "Model_comments" Description: ""
--     * Slot: Model_id Description: Autocreated FK slot
--     * Slot: comments Description: Comments about the model
-- # Class: "EvidenceItem_with_objects" Description: ""
--     * Slot: EvidenceItem_id Description: Autocreated FK slot
--     * Slot: with_objects_id Description: Supporting database entities or terms
-- # Class: "EnabledByProteinComplexAssociation_members" Description: ""
--     * Slot: EnabledByProteinComplexAssociation_id Description: Autocreated FK slot
--     * Slot: members_id Description: The gene products that are part of the complex

CREATE TABLE "Activity" (
	id TEXT NOT NULL, 
	"Model_id" TEXT, 
	enabled_by_id INTEGER, 
	molecular_function_id INTEGER, 
	occurs_in_id INTEGER, 
	part_of_id INTEGER, 
	has_primary_input_id INTEGER, 
	has_primary_output_id INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("Model_id") REFERENCES "Model" (id), 
	FOREIGN KEY(enabled_by_id) REFERENCES "EnabledByAssociation" (id), 
	FOREIGN KEY(molecular_function_id) REFERENCES "MolecularFunctionAssociation" (id), 
	FOREIGN KEY(occurs_in_id) REFERENCES "CellularAnatomicalEntityAssociation" (id), 
	FOREIGN KEY(part_of_id) REFERENCES "BiologicalProcessAssociation" (id), 
	FOREIGN KEY(has_primary_input_id) REFERENCES "MoleculeAssociation" (id), 
	FOREIGN KEY(has_primary_output_id) REFERENCES "MoleculeAssociation" (id)
);
CREATE TABLE "Association" (
	id INTEGER NOT NULL, 
	type TEXT, 
	PRIMARY KEY (id)
);
CREATE TABLE "MoleculeAssociation" (
	id INTEGER NOT NULL, 
	term TEXT, 
	type TEXT, 
	"Activity_id" TEXT, 
	PRIMARY KEY (id), 
	FOREIGN KEY(term) REFERENCES "MoleculeTermObject" (id), 
	FOREIGN KEY("Activity_id") REFERENCES "Activity" (id)
);
CREATE TABLE "TermObject" (
	id TEXT NOT NULL, 
	label TEXT, 
	type TEXT, 
	obsolete BOOLEAN, 
	PRIMARY KEY (id)
);
CREATE TABLE "PublicationObject" (
	abstract_text TEXT, 
	full_text TEXT, 
	id TEXT NOT NULL, 
	label TEXT, 
	type TEXT, 
	obsolete BOOLEAN, 
	PRIMARY KEY (id)
);
CREATE TABLE "EvidenceTermObject" (
	id TEXT NOT NULL, 
	label TEXT, 
	type TEXT, 
	obsolete BOOLEAN, 
	PRIMARY KEY (id)
);
CREATE TABLE "MolecularFunctionTermObject" (
	id TEXT NOT NULL, 
	label TEXT, 
	type TEXT, 
	obsolete BOOLEAN, 
	PRIMARY KEY (id)
);
CREATE TABLE "BiologicalProcessTermObject" (
	id TEXT NOT NULL, 
	label TEXT, 
	type TEXT, 
	obsolete BOOLEAN, 
	PRIMARY KEY (id)
);
CREATE TABLE "CellularAnatomicalEntityTermObject" (
	id TEXT NOT NULL, 
	label TEXT, 
	type TEXT, 
	obsolete BOOLEAN, 
	PRIMARY KEY (id)
);
CREATE TABLE "MoleculeTermObject" (
	id TEXT NOT NULL, 
	label TEXT, 
	type TEXT, 
	obsolete BOOLEAN, 
	PRIMARY KEY (id)
);
CREATE TABLE "CellTypeTermObject" (
	id TEXT NOT NULL, 
	label TEXT, 
	type TEXT, 
	obsolete BOOLEAN, 
	PRIMARY KEY (id)
);
CREATE TABLE "GrossAnatomicalStructureTermObject" (
	id TEXT NOT NULL, 
	label TEXT, 
	type TEXT, 
	obsolete BOOLEAN, 
	PRIMARY KEY (id)
);
CREATE TABLE "PhaseTermObject" (
	id TEXT NOT NULL, 
	label TEXT, 
	type TEXT, 
	obsolete BOOLEAN, 
	PRIMARY KEY (id)
);
CREATE TABLE "InformationBiomacromoleculeTermObject" (
	id TEXT NOT NULL, 
	label TEXT, 
	type TEXT, 
	obsolete BOOLEAN, 
	PRIMARY KEY (id)
);
CREATE TABLE "GeneProductTermObject" (
	id TEXT NOT NULL, 
	label TEXT, 
	type TEXT, 
	obsolete BOOLEAN, 
	PRIMARY KEY (id)
);
CREATE TABLE "ProteinComplexTermObject" (
	id TEXT NOT NULL, 
	label TEXT, 
	type TEXT, 
	obsolete BOOLEAN, 
	PRIMARY KEY (id)
);
CREATE TABLE "TaxonTermObject" (
	id TEXT NOT NULL, 
	label TEXT, 
	type TEXT, 
	obsolete BOOLEAN, 
	PRIMARY KEY (id)
);
CREATE TABLE "PredicateTermObject" (
	id TEXT NOT NULL, 
	label TEXT, 
	type TEXT, 
	obsolete BOOLEAN, 
	PRIMARY KEY (id)
);
CREATE TABLE "Model" (
	id TEXT NOT NULL, 
	title TEXT, 
	taxon TEXT, 
	status VARCHAR(11), 
	PRIMARY KEY (id), 
	FOREIGN KEY(taxon) REFERENCES "TaxonTermObject" (id)
);
CREATE TABLE "EnabledByAssociation" (
	id INTEGER NOT NULL, 
	term TEXT, 
	type TEXT, 
	PRIMARY KEY (id), 
	FOREIGN KEY(term) REFERENCES "InformationBiomacromoleculeTermObject" (id)
);
CREATE TABLE "EnabledByGeneProductAssociation" (
	id INTEGER NOT NULL, 
	term TEXT, 
	type TEXT, 
	PRIMARY KEY (id), 
	FOREIGN KEY(term) REFERENCES "GeneProductTermObject" (id)
);
CREATE TABLE "EnabledByProteinComplexAssociation" (
	id INTEGER NOT NULL, 
	term TEXT, 
	type TEXT, 
	PRIMARY KEY (id), 
	FOREIGN KEY(term) REFERENCES "ProteinComplexTermObject" (id)
);
CREATE TABLE "CausalAssociation" (
	id INTEGER NOT NULL, 
	predicate TEXT, 
	downstream_activity TEXT, 
	type TEXT, 
	"Activity_id" TEXT, 
	PRIMARY KEY (id), 
	FOREIGN KEY(predicate) REFERENCES "PredicateTermObject" (id), 
	FOREIGN KEY(downstream_activity) REFERENCES "Activity" (id), 
	FOREIGN KEY("Activity_id") REFERENCES "Activity" (id)
);
CREATE TABLE "TermAssociation" (
	id INTEGER NOT NULL, 
	term TEXT, 
	type TEXT, 
	PRIMARY KEY (id), 
	FOREIGN KEY(term) REFERENCES "TermObject" (id)
);
CREATE TABLE "MolecularFunctionAssociation" (
	id INTEGER NOT NULL, 
	term TEXT, 
	type TEXT, 
	PRIMARY KEY (id), 
	FOREIGN KEY(term) REFERENCES "MolecularFunctionTermObject" (id)
);
CREATE TABLE "BiologicalProcessAssociation" (
	id INTEGER NOT NULL, 
	happens_during TEXT, 
	part_of TEXT, 
	term TEXT, 
	type TEXT, 
	PRIMARY KEY (id), 
	FOREIGN KEY(happens_during) REFERENCES "PhaseTermObject" (id), 
	FOREIGN KEY(part_of) REFERENCES "BiologicalProcessTermObject" (id), 
	FOREIGN KEY(term) REFERENCES "BiologicalProcessTermObject" (id)
);
CREATE TABLE "GrossAnatomyAssociation" (
	id INTEGER NOT NULL, 
	term TEXT, 
	type TEXT, 
	part_of_id INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY(term) REFERENCES "GrossAnatomicalStructureTermObject" (id), 
	FOREIGN KEY(part_of_id) REFERENCES "GrossAnatomyAssociation" (id)
);
CREATE TABLE "CellTypeAssociation" (
	id INTEGER NOT NULL, 
	term TEXT, 
	type TEXT, 
	part_of_id INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY(term) REFERENCES "CellTypeTermObject" (id), 
	FOREIGN KEY(part_of_id) REFERENCES "GrossAnatomyAssociation" (id)
);
CREATE TABLE "Object" (
	id TEXT NOT NULL, 
	label TEXT, 
	type TEXT, 
	obsolete BOOLEAN, 
	"Model_id" TEXT, 
	PRIMARY KEY (id), 
	FOREIGN KEY("Model_id") REFERENCES "Model" (id)
);
CREATE TABLE "Model_comments" (
	"Model_id" TEXT, 
	comments TEXT, 
	PRIMARY KEY ("Model_id", comments), 
	FOREIGN KEY("Model_id") REFERENCES "Model" (id)
);
CREATE TABLE "EnabledByProteinComplexAssociation_members" (
	"EnabledByProteinComplexAssociation_id" INTEGER, 
	members_id TEXT, 
	PRIMARY KEY ("EnabledByProteinComplexAssociation_id", members_id), 
	FOREIGN KEY("EnabledByProteinComplexAssociation_id") REFERENCES "EnabledByProteinComplexAssociation" (id), 
	FOREIGN KEY(members_id) REFERENCES "GeneProductTermObject" (id)
);
CREATE TABLE "CellularAnatomicalEntityAssociation" (
	id INTEGER NOT NULL, 
	term TEXT, 
	type TEXT, 
	part_of_id INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY(term) REFERENCES "CellularAnatomicalEntityTermObject" (id), 
	FOREIGN KEY(part_of_id) REFERENCES "CellTypeAssociation" (id)
);
CREATE TABLE "EvidenceItem" (
	id INTEGER NOT NULL, 
	term TEXT, 
	reference TEXT, 
	"Association_id" INTEGER, 
	"EnabledByAssociation_id" INTEGER, 
	"EnabledByGeneProductAssociation_id" INTEGER, 
	"EnabledByProteinComplexAssociation_id" INTEGER, 
	"CausalAssociation_id" INTEGER, 
	"TermAssociation_id" INTEGER, 
	"MolecularFunctionAssociation_id" INTEGER, 
	"BiologicalProcessAssociation_id" INTEGER, 
	"CellularAnatomicalEntityAssociation_id" INTEGER, 
	"CellTypeAssociation_id" INTEGER, 
	"GrossAnatomyAssociation_id" INTEGER, 
	"MoleculeAssociation_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY(term) REFERENCES "EvidenceTermObject" (id), 
	FOREIGN KEY(reference) REFERENCES "PublicationObject" (id), 
	FOREIGN KEY("Association_id") REFERENCES "Association" (id), 
	FOREIGN KEY("EnabledByAssociation_id") REFERENCES "EnabledByAssociation" (id), 
	FOREIGN KEY("EnabledByGeneProductAssociation_id") REFERENCES "EnabledByGeneProductAssociation" (id), 
	FOREIGN KEY("EnabledByProteinComplexAssociation_id") REFERENCES "EnabledByProteinComplexAssociation" (id), 
	FOREIGN KEY("CausalAssociation_id") REFERENCES "CausalAssociation" (id), 
	FOREIGN KEY("TermAssociation_id") REFERENCES "TermAssociation" (id), 
	FOREIGN KEY("MolecularFunctionAssociation_id") REFERENCES "MolecularFunctionAssociation" (id), 
	FOREIGN KEY("BiologicalProcessAssociation_id") REFERENCES "BiologicalProcessAssociation" (id), 
	FOREIGN KEY("CellularAnatomicalEntityAssociation_id") REFERENCES "CellularAnatomicalEntityAssociation" (id), 
	FOREIGN KEY("CellTypeAssociation_id") REFERENCES "CellTypeAssociation" (id), 
	FOREIGN KEY("GrossAnatomyAssociation_id") REFERENCES "GrossAnatomyAssociation" (id), 
	FOREIGN KEY("MoleculeAssociation_id") REFERENCES "MoleculeAssociation" (id)
);
CREATE TABLE "ProvenanceInfo" (
	id INTEGER NOT NULL, 
	contributor TEXT, 
	created TEXT, 
	date TEXT, 
	provided_by TEXT, 
	"Model_id" TEXT, 
	"Activity_id" TEXT, 
	"EvidenceItem_id" INTEGER, 
	"Association_id" INTEGER, 
	"EnabledByAssociation_id" INTEGER, 
	"EnabledByGeneProductAssociation_id" INTEGER, 
	"EnabledByProteinComplexAssociation_id" INTEGER, 
	"CausalAssociation_id" INTEGER, 
	"TermAssociation_id" INTEGER, 
	"MolecularFunctionAssociation_id" INTEGER, 
	"BiologicalProcessAssociation_id" INTEGER, 
	"CellularAnatomicalEntityAssociation_id" INTEGER, 
	"CellTypeAssociation_id" INTEGER, 
	"GrossAnatomyAssociation_id" INTEGER, 
	"MoleculeAssociation_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("Model_id") REFERENCES "Model" (id), 
	FOREIGN KEY("Activity_id") REFERENCES "Activity" (id), 
	FOREIGN KEY("EvidenceItem_id") REFERENCES "EvidenceItem" (id), 
	FOREIGN KEY("Association_id") REFERENCES "Association" (id), 
	FOREIGN KEY("EnabledByAssociation_id") REFERENCES "EnabledByAssociation" (id), 
	FOREIGN KEY("EnabledByGeneProductAssociation_id") REFERENCES "EnabledByGeneProductAssociation" (id), 
	FOREIGN KEY("EnabledByProteinComplexAssociation_id") REFERENCES "EnabledByProteinComplexAssociation" (id), 
	FOREIGN KEY("CausalAssociation_id") REFERENCES "CausalAssociation" (id), 
	FOREIGN KEY("TermAssociation_id") REFERENCES "TermAssociation" (id), 
	FOREIGN KEY("MolecularFunctionAssociation_id") REFERENCES "MolecularFunctionAssociation" (id), 
	FOREIGN KEY("BiologicalProcessAssociation_id") REFERENCES "BiologicalProcessAssociation" (id), 
	FOREIGN KEY("CellularAnatomicalEntityAssociation_id") REFERENCES "CellularAnatomicalEntityAssociation" (id), 
	FOREIGN KEY("CellTypeAssociation_id") REFERENCES "CellTypeAssociation" (id), 
	FOREIGN KEY("GrossAnatomyAssociation_id") REFERENCES "GrossAnatomyAssociation" (id), 
	FOREIGN KEY("MoleculeAssociation_id") REFERENCES "MoleculeAssociation" (id)
);
CREATE TABLE "EvidenceItem_with_objects" (
	"EvidenceItem_id" INTEGER, 
	with_objects_id TEXT, 
	PRIMARY KEY ("EvidenceItem_id", with_objects_id), 
	FOREIGN KEY("EvidenceItem_id") REFERENCES "EvidenceItem" (id), 
	FOREIGN KEY(with_objects_id) REFERENCES "Object" (id)
);