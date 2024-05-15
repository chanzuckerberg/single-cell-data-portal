DIMENSION_NAME_MAP_CENSUS_TO_WMG = {
    "tissue_ontology_term_id": "tissue_original_ontology_term_id",
    "tissue_general_ontology_term_id": "tissue_ontology_term_id",
}

EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG = "expression_summary_and_cell_counts_cube_created"
EXPRESSION_SUMMARY_DIFFEXP_CUBES_CREATED_FLAG = "expression_summary_diffexp_cubes_created"
EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG = "expression_summary_default_cube_created"
MARKER_GENES_CUBE_CREATED_FLAG = "marker_genes_cube_created"
FILTER_RELATIONSHIPS_CREATED_FLAG = "filter_relationships_created"
PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG = "primary_filter_dimensions_created"
DATASET_METADATA_CREATED_FLAG = "dataset_metadata_created"
CELL_TYPE_ANCESTORS_CREATED_FLAG = "cell_type_ancestors_created"
CELL_TYPE_ORDERING_CREATED_FLAG = "cell_type_ordering_created"

PIPELINE_STATE_FILENAME = "pipeline_state.json"

# Minimum number of expressed genes for a cell to be included in the corpus.
# See the following document for further details:
# https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/scExpression/scExpression-documentation.md#removal-of-low-coverage-cells
GENE_EXPRESSION_COUNT_MIN_THRESHOLD = 500

# Minimum value for raw expression counts that will be used to filter out computed RankIt values. Details:
# https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/scExpression/scExpression-documentation.md#removal-of-noisy-ultra-low-expression-values
NORM_EXPR_COUNT_FILTERING_MIN_THRESHOLD = 3

ASSAYS_FOR_GENE_LENGTH_NORMALIZATION = [
    "EFO:0008930",  # Smart-seq
    "EFO:0008931",  # Smart-seq2
    "EFO:0700016",  # Smart-seq v4
]

TARGET_LIBRARY_SIZE = 10_000

WMG_DATA_SCHEMA_VERSION = "v4"

MAXIMUM_ADMISSIBLE_CENSUS_SCHEMA_MAJOR_VERSION = 2

HIGH_LEVEL_TISSUES = [
    "UBERON:0000178",  # blood
    "UBERON:0002048",  # lung
    "UBERON:0002106",  # spleen
    "UBERON:0002371",  # bone marrow
    "UBERON:0002107",  # liver
    "UBERON:0002113",  # kidney
    "UBERON:0000955",  # brain
    "UBERON:0002240",  # spinal cord
    "UBERON:0000310",  # breast
    "UBERON:0000948",  # heart
    "UBERON:0002097",  # skin of body
    "UBERON:0000970",  # eye
    "UBERON:0001264",  # pancreas
    "UBERON:0001043",  # esophagus
    "UBERON:0001155",  # colon
    "UBERON:0000059",  # large intestine
    "UBERON:0002108",  # small intestine
    "UBERON:0000160",  # intestine
    "UBERON:0000945",  # stomach
    "UBERON:0001836",  # saliva
    "UBERON:0001723",  # tongue
    "UBERON:0001013",  # adipose tissue
    "UBERON:0000473",  # testis
    "UBERON:0002367",  # prostate gland
    "UBERON:0000057",  # urethra
    "UBERON:0000056",  # ureter
    "UBERON:0003889",  # fallopian tube
    "UBERON:0000995",  # uterus
    "UBERON:0000992",  # ovary
    "UBERON:0002110",  # gall bladder
    "UBERON:0001255",  # urinary bladder
    "UBERON:0018707",  # bladder organ
    "UBERON:0000922",  # embryo
    "UBERON:0004023",  # ganglionic eminence --> this a part of the embryo, remove in case generality is desired
    "UBERON:0001987",  # placenta
    "UBERON:0007106",  # chorionic villus
    "UBERON:0002369",  # adrenal gland
    "UBERON:0002368",  # endocrine gland
    "UBERON:0002365",  # exocrine gland
    "UBERON:0000030",  # lamina propria
    "UBERON:0000029",  # lymph node
    "UBERON:0004536",  # lymph vasculature
    "UBERON:0001015",  # musculature
    "UBERON:0000004",  # nose
    "UBERON:0003688",  # omentum
    "UBERON:0000977",  # pleura
    "UBERON:0002370",  # thymus
    "UBERON:0002049",  # vasculature
    "UBERON:0009472",  # axilla
    "UBERON:0001087",  # pleural fluid
    "UBERON:0000344",  # mucosa
    "UBERON:0001434",  # skeletal system
    "UBERON:0002228",  # rib
    "UBERON:0003129",  # skull
    "UBERON:0004537",  # blood vasculature
    "UBERON:0002405",  # immune system
    "UBERON:0001009",  # circulatory system
    "UBERON:0001007",  # digestive system
    "UBERON:0001017",  # central nervous system
    "UBERON:0001008",  # renal system
    "UBERON:0000990",  # reproductive system
    "UBERON:0001004",  # respiratory system
    "UBERON:0000010",  # peripheral nervous system
    "UBERON:0001032",  # sensory system
    "UBERON:0002046",  # thyroid gland
    "UBERON:0004535",  # cardiovascular system
    "UBERON:0000949",  # endocrine system
    "UBERON:0002330",  # exocrine system
    "UBERON:0002390",  # hematopoietic system
    "UBERON:0000383",  # musculature of body
    "UBERON:0001465",  # knee
    "UBERON:0001016",  # nervous system
    "UBERON:0001348",  # brown adipose tissue
    "UBERON:0015143",  # mesenteric fat pad
    "UBERON:0000175",  # pleural effusion
    "UBERON:0001416",  # skin of abdomen
    "UBERON:0001868",  # skin of chest
    "UBERON:0001511",  # skin of leg
    "UBERON:0002190",  # subcutaneous adipose tissue
    "UBERON:0000014",  # zone of skin
    "UBERON:0000916",  # abdomen
]
