# Derived from:
# https://github.com/geneontology/web-components/blob/5d87e593121eafe6ac4690fa4591f88aa5a03fd8/packages/web-components/src/globals/%40noctua.form/data/taxon-dataset.json
SPECIES_CODES = [
    "Atal"
    "Btau"
    "Cele"
    "Cfam"
    "Ddis"
    "Dmel"
    "Drer"
    "Ggal"
    "Hsap"
    "Mmus"
    "Pseudomonas"
    "Rnor"
    "Scer"
    "Sjap"
    "Solanaceae"
    "Spom"
    "Sscr"
    "Xenopus"
]


def remove_species_code_suffix(label: str) -> str:
    for code in SPECIES_CODES:
        label = label.removesuffix(code).strip()
    return label
