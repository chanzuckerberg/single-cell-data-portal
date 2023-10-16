DIMENSION_NAME_MAP_CENSUS_TO_WMG = {
    "tissue_ontology_term_id": "tissue_original_ontology_term_id",
    "tissue_general_ontology_term_id": "tissue_ontology_term_id",
}

EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG = "expression_summary_and_cell_counts_cube_created"
EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG = "expression_summary_default_cube_created"
MARKER_GENES_CUBE_CREATED_FLAG = "marker_genes_cube_created"
FILTER_RELATIONSHIPS_CREATED_FLAG = "filter_relationships_created"
PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG = "primary_filter_dimensions_created"
DATASET_METADATA_CREATED_FLAG = "dataset_metadata_created"

PIPELINE_STATE_FILENAME = "pipeline_state.json"

# Minimum number of expressed genes for a cell to be included in the corpus.
# See the following document for further details:
# https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/scExpression/scExpression-documentation.md#removal-of-low-coverage-cells
GENE_EXPRESSION_COUNT_MIN_THRESHOLD = 500

# Minimum value for raw expression counts that will be used to filter out computed RankIt values. Details:
# https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/scExpression/scExpression-documentation.md#removal-of-noisy-ultra-low-expression-values
RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD = 2
