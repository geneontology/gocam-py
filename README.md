# gocam

GO CAM Data Model (Python)

## Website

[https://geneontology.github.io/gocam-py](https://geneontology.github.io/gocam-py)

## About

Example yaml/json

```yaml
id: gomodel:568b0f9600000284
title: Antibacterial innate immune response in the intestine via MAPK cascade (C.
  elegans)
taxon: NCBITaxon:6239
status: production
comments:
- 'Automated change 2023-03-16: RO:0002212 replaced by RO:0002630'
- 'Automated change 2023-03-16: RO:0002213 replaced by RO:0002629'
activities:
- id: gomodel:568b0f9600000284/57ec3a7e00000079
  enabled_by: WB:WBGene00006575
  molecular_function:
    evidence:
    - term: ECO:0000314
      reference: PMID:15625192
      provenances:
      - contributor: https://orcid.org/0000-0002-1706-4196
        date: '2019-09-23'
    provenances: []
    term: GO:0035591
  occurs_in:
    evidence:
    - term: ECO:0000314
      reference: PMID:15625192
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2019-06-28'
    term: GO:0005737
  part_of:
    evidence:
    - term: ECO:0000315
      reference: PMID:19837372
      with_objects:
      - WB:WBVar00241222|WB:WBVar00241223
      provenances:
      - contributor: https://orcid.org/0000-0002-1706-4196
        date: '2021-07-08'
    term: GO:0140367
  causal_associations:
  - evidence:
    - term: ECO:0000315
      reference: PMID:15123841
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2019-05-31'
    predicate: RO:0002629
    downstream_activity: gomodel:568b0f9600000284/57ec3a7e00000109
- id: gomodel:568b0f9600000284/580685bd00000135
  enabled_by: WB:WBGene00011979
  molecular_function:
    evidence:
    - term: ECO:0000307
      reference: GO_REF:0000015
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-11-12'
    provenances: []
    term: GO:0003674
  occurs_in:
    evidence:
    - term: ECO:0000307
      reference: GO_REF:0000015
      provenances:
      - contributor: https://orcid.org/0000-0002-1706-4196
        date: '2021-11-29'
    term: GO:0110165
  part_of:
    evidence:
    - term: ECO:0000270
      reference: PMID:17096597
      provenances:
      - contributor: https://orcid.org/0000-0002-1706-4196
        date: '2021-07-08'
    term: GO:0140367
- id: gomodel:568b0f9600000284/5745387b00000588
  enabled_by: WB:WBGene00012019
  molecular_function:
    evidence:
    - term: ECO:0000314
      reference: PMID:17728253
      provenances:
      - contributor: https://orcid.org/0000-0002-1706-4196
        date: '2019-09-23'
    provenances: []
    term: GO:0004674
  occurs_in:
    evidence:
    - term: ECO:0000314
      reference: PMID:17728253
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2019-07-01'
    term: GO:0009898
  part_of:
    evidence:
    - term: ECO:0000315
      reference: PMID:19371715
      provenances:
      - contributor: https://orcid.org/0000-0002-1706-4196
        date: '2016-10-19'
    term: GO:0140367
  causal_associations:
  - evidence:
    - term: ECO:0000315
      reference: PMID:19371715
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2019-07-02'
    predicate: RO:0002629
    downstream_activity: gomodel:568b0f9600000284/568b0f9600000285
- id: gomodel:568b0f9600000284/5b528b1100002286
  enabled_by: WB:WBGene00006599
  molecular_function:
    evidence:
    - term: ECO:0000501
      reference: GO_REF:0000037
      with_objects:
      - UniProtKB-KW:KW-0723
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-09-06'
    - term: ECO:0000315
      reference: PMID:23072806
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-09-06'
    - term: ECO:0000501
      reference: GO_REF:0000002
      with_objects:
      - InterPro:IPR000961
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-09-06'
    - term: ECO:0000318
      reference: PAINT_REF:24356
      with_objects:
      - PANTHER:PTN000683254
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-09-06'
    provenances: []
    term: GO:0004674
  part_of:
    evidence:
    - term: ECO:0000315
      reference: PMID:22470487
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-09-06'
    term: GO:0002225
  causal_associations:
  - evidence:
    - term: ECO:0000316
      reference: PMID:19371715
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2019-07-01'
    predicate: RO:0002629
    downstream_activity: gomodel:568b0f9600000284/5745387b00000588
- id: gomodel:568b0f9600000284/57ec3a7e00000119
  enabled_by: WB:WBGene00004758
  molecular_function:
    evidence:
    - term: ECO:0000314
      reference: PMID:11751572
      provenances:
      - contributor: https://orcid.org/0000-0002-1706-4196
        date: '2019-09-23'
    provenances: []
    term: GO:0004708
  occurs_in:
    evidence:
    - term: ECO:0000318
      reference: PMID:21873635
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2019-06-28'
    term: GO:0005737
  part_of:
    evidence:
    - term: ECO:0000315
      reference: PMID:12142542
      provenances:
      - contributor: https://orcid.org/0000-0002-1706-4196
        date: '2020-10-28'
    term: GO:0140367
  causal_associations:
  - evidence:
    - term: ECO:0000315
      reference: PMID:12142542
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2019-05-31'
    predicate: RO:0002629
    downstream_activity: gomodel:568b0f9600000284/568b0f9600000285
- id: gomodel:568b0f9600000284/568b0f9600000285
  enabled_by: WB:WBGene00004055
  molecular_function:
    evidence:
    - term: ECO:0000314
      reference: PMID:20369020
      provenances:
      - contributor: https://orcid.org/0000-0002-1706-4196
        date: '2020-10-28'
    provenances: []
    term: GO:0004707
  occurs_in:
    evidence:
    - term: ECO:0000314
      reference: PMID:20133945
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2019-07-01'
    term: GO:0005829
  part_of:
    evidence:
    - term: ECO:0000315
      reference: PMID:12142542
      provenances:
      - contributor: https://orcid.org/0000-0002-1706-4196
        date: '2020-10-28'
    term: GO:0140367
  causal_associations:
  - evidence:
    - term: ECO:0000314
      reference: PMID:20369020
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2019-05-31'
    predicate: RO:0002629
    downstream_activity: gomodel:568b0f9600000284/568b0f9600000287
- id: gomodel:568b0f9600000284/57ec3a7e00000109
  enabled_by: WB:WBGene00003822
  molecular_function:
    evidence:
    - term: ECO:0000314
      reference: PMID:11751572
      provenances:
      - contributor: https://orcid.org/0000-0002-1706-4196
        date: '2019-09-23'
    provenances: []
    term: GO:0004709
  occurs_in:
    evidence:
    - term: ECO:0000314
      reference: PMID:15625192
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2019-06-28'
    term: GO:0005737
  part_of:
    evidence:
    - term: ECO:0000315
      reference: PMID:12142542
      provenances:
      - contributor: https://orcid.org/0000-0002-1706-4196
        date: '2020-10-28'
    term: GO:0140367
  causal_associations:
  - evidence:
    - term: ECO:0000315
      reference: PMID:12142542
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2019-05-31'
    predicate: RO:0002629
    downstream_activity: gomodel:568b0f9600000284/57ec3a7e00000119
- id: gomodel:568b0f9600000284/5b91dbd100000652
  enabled_by: WB:WBGene00006923
  molecular_function:
    evidence:
    - term: ECO:0000314
      reference: PMID:15116070
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-10-08'
    provenances: []
    term: GO:0017017
  occurs_in:
    evidence:
    - term: ECO:0000250
      reference: PMID:15116070
      with_objects:
      - UniProt:Q920R2
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2019-07-01'
    term: GO:0005737
  part_of:
    evidence:
    - term: ECO:0000316
      reference: PMID:22554143
      with_objects:
      - WB:WBGene00002187
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-10-08'
    - term: ECO:0000315
      reference: PMID:22554143
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-10-08'
    term: GO:1900425
  causal_associations:
  - evidence:
    - term: ECO:0000315
      reference: PMID:22554143
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-10-09'
    predicate: RO:0002630
    downstream_activity: gomodel:568b0f9600000284/5b91dbd100000659
  - evidence:
    - term: ECO:0000315
      reference: PMID:15256594
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2019-06-20'
    - term: ECO:0000315
      reference: PMID:22554143
      provenances:
      - contributor: https://orcid.org/0000-0002-1706-4196
        date: '2020-10-28'
    predicate: RO:0002630
    downstream_activity: gomodel:568b0f9600000284/568b0f9600000285
- id: gomodel:568b0f9600000284/5b91dbd100000659
  enabled_by: WB:WBGene00002187
  molecular_function:
    evidence:
    - term: ECO:0000353
      reference: PMID:12435362
      with_objects:
      - WB:WBGene00001599
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-10-08'
    - term: ECO:0000353
      reference: PMID:12435362
      with_objects:
      - WB:WBGene00001600
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-10-08'
    - term: ECO:0000353
      reference: PMID:23437011
      with_objects:
      - WB:WBGene00001345
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-10-08'
    - term: ECO:0000353
      reference: PMID:23437011
      with_objects:
      - WB:WBGene00001345
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-10-08'
    provenances: []
    term: GO:0005515
  occurs_in:
    evidence:
    - term: ECO:0000314
      reference: PMID:17699606
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2019-07-01'
    term: GO:0005737
  part_of:
    evidence:
    - term: ECO:0000315
      reference: PMID:22554143
      with_objects:
      - WB:WBGene00002187
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-10-08'
    term: GO:1900181
- id: gomodel:568b0f9600000284/568b0f9600000287
  enabled_by: WB:WBGene00000223
  molecular_function:
    evidence:
    - term: ECO:0000250
      reference: PMID:20369020
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-10-29'
    provenances: []
    term: GO:0000981
  occurs_in:
    evidence:
    - term: ECO:0000314
      reference: PMID:20369020
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-10-29'
    term: GO:0005634
  part_of:
    evidence:
    - term: ECO:0000315
      reference: PMID:20369020
      provenances:
      - contributor: https://orcid.org/0000-0002-1706-4196
        date: '2020-10-28'
    term: GO:0140367
  causal_associations:
  - evidence:
    - term: ECO:0000315
      reference: PMID:20369020
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-09-26'
    - term: ECO:0000315
      reference: PMID:19837372
      provenances:
      - contributor: https://orcid.org/0000-0002-3013-9906
        date: '2018-10-31'
    predicate: RO:0002629
    downstream_activity: gomodel:568b0f9600000284/580685bd00000135
objects:
- id: WB:WBGene00006599
  label: tpa-1 Cele
- id: GO:0004674
  label: protein serine/threonine kinase activity
- id: GO:0002225
  label: positive regulation of antimicrobial peptide production
- id: ECO:0000501
  label: evidence used in automatic assertion
- id: ECO:0000315
  label: mutant phenotype evidence used in manual assertion
- id: ECO:0000318
  label: biological aspect of ancestor evidence used in manual assertion
- id: WB:WBGene00006923
  label: vhp-1 Cele
- id: GO:0017017
  label: MAP kinase tyrosine/serine/threonine phosphatase activity
- id: GO:1900425
  label: negative regulation of defense response to bacterium
- id: ECO:0000314
  label: direct assay evidence used in manual assertion
- id: ECO:0000316
  label: genetic interaction evidence used in manual assertion
- id: WB:WBGene00002187
  label: kgb-1 Cele
- id: GO:0005515
  label: protein binding
- id: GO:1900181
  label: negative regulation of protein localization to nucleus
- id: ECO:0000353
  label: physical interaction evidence used in manual assertion
- id: ECO:0000250
  label: sequence similarity evidence used in manual assertion
- id: ECO:0000307
  label: no evidence data found used in manual assertion
- id: GO:0005737
  label: cytoplasm
- id: GO:0005829
  label: cytosol
- id: GO:0009898
  label: cytoplasmic side of plasma membrane
- id: WB:WBGene00000223
  label: atf-7 Cele
- id: WB:WBGene00004055
  label: pmk-1 Cele
- id: WB:WBGene00004758
  label: sek-1 Cele
- id: GO:0004707
  label: MAP kinase activity
- id: GO:0000981
  label: DNA-binding transcription factor activity, RNA polymerase II-specific
- id: GO:0005634
  label: nucleus
- id: GO:0016045
  label: detection of bacterium
- id: GO:0003674
  label: molecular_function
- id: GO:0035591
  label: signaling adaptor activity
- id: WB:WBGene00006575
  label: tir-1 Cele
- id: GO:0004709
  label: MAP kinase kinase kinase activity
- id: WB:WBGene00003822
  label: nsy-1 Cele
- id: GO:0004708
  label: MAP kinase kinase activity
- id: WB:WBGene00012019
  label: dkf-2 Cele
- id: GO:0140367
  label: antibacterial innate immune response
- id: WB:WBGene00011979
  label: sysm-1 Cele
- id: ECO:0000270
  label: expression pattern evidence used in manual assertion
- id: GO:0110165
  label: cellular anatomical entity

```

