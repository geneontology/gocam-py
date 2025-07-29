"""
Taxonomy vocabularies and utilities for GOCAM
"""

class TaxonVocabulary:
    """
    Vocabulary of taxonomy terms and constants used in GOCAM
    """
    # Host/model organism taxa that should be prioritized as primary
    HOST_TAXA = [
        "NCBITaxon:9606",  # Human
        "NCBITaxon:10090", # Mouse
        "NCBITaxon:10116", # Rat
        "NCBITaxon:7227",  # Fruit fly
        "NCBITaxon:6239",  # C. elegans
        "NCBITaxon:7955",  # Zebrafish
        "NCBITaxon:4932",  # S. cerevisiae
        "NCBITaxon:3702",  # Arabidopsis
    ]
    
    # Annotation key for taxa in models
    TAXON_ANNOTATION_KEY = "in_taxon"  # This is the normalized version used in the code (after _normalize_property)
    BIOLINK_TAXON_KEY = "https://w3id.org/biolink/vocab/in_taxon"  # Original biolink key
    LEGACY_TAXON_KEY = "in_taxon"     # Same as normalized for now
    
    @classmethod
    def is_host_taxon(cls, taxon_id):
        """
        Check if a taxon is a host/model organism
        
        :param taxon_id: Taxon ID in CURIE format (e.g., "NCBITaxon:9606")
        :return: True if the taxon is a host organism, False otherwise
        """
        return taxon_id in cls.HOST_TAXA