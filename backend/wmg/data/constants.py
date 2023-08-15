# Minimum number of expressed genes for a cell to be included in the corpus.
# See the following document for further details:
# https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/scExpression/scExpression-documentation.md#removal-of-low-coverage-cells
GENE_EXPRESSION_COUNT_MIN_THRESHOLD = 500

# Minimum value for raw expression counts that will be used to filter out computed RankIt values. Details:
# https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/scExpression/scExpression-documentation.md#removal-of-noisy-ultra-low-expression-values
RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD = 2

# this mirrors the assays used by the census builder
# https://github.com/chanzuckerberg/cellxgene-census/blob/main/tools/cellxgene_census_builder/src/cellxgene_census_builder/build_soma/globals.py
INCLUDED_ASSAYS = {
    "EFO:0008720": "DroNc-seq",
    "EFO:0008722": "Drop-seq",
    "EFO:0008780": "inDrop",
    "EFO:0008913": "single-cell RNA sequencing",
    "EFO:0008919": "Seq-Well",
    "EFO:0008930": "Smart-seq",
    "EFO:0008931": "Smart-seq2",
    "EFO:0008953": "STRT-seq",
    "EFO:0008995": "10x technology",
    "EFO:0009899": "10x 3' v2",
    "EFO:0009900": "10x 5' v2",
    "EFO:0009901": "10x 3' v1",
    "EFO:0009922": "10x 3' v3",
    "EFO:0010010": "CEL-seq2",
    "EFO:0010183": "single cell library construction",
    "EFO:0010550": "sci-RNA-seq",
    "EFO:0011025": "10x 5' v1",
    "EFO:0030002": "microwell-seq",
    "EFO:0030003": "10x 3' transcription profiling",
    "EFO:0030004": "10x 5' transcription profiling",
    "EFO:0030019": "Seq-Well S3",
    "EFO:0700003": "BD Rhapsody Whole Transcriptome Analysis",
    "EFO:0700004": "BD Rhapsody Targeted mRNA",
}


CL_BASIC_PERMANENT_URL_PRONTO = "https://github.com/obophenotype/cell-ontology/releases/latest/download/cl-basic.obo"
CL_BASIC_PERMANENT_URL_OWL = "https://github.com/obophenotype/cell-ontology/releases/latest/download/cl-basic.owl"
