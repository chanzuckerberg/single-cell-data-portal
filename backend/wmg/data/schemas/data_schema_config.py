WMG_DATA_SCHEMA_VERSION = "v1"

# These are the valid attributes and dimensions consulted by the
# wmg pipeline (writer) to determine the list of attributes and dimensions
# TO RETRIEVE from the cube when performing a query.
# This is important because the cube allocates a large amount of memory for
# each attribute and dimension retrieved. Therefore, strictly specifying the
# attributes and dimensions to retrieve makes the query more memory efficient
WRITER_WMG_CUBE_QUERY_VALID_ATTRIBUTES = ["gene_ontology_term_id", "cell_type_ontology_term_id"]
WRITER_WMG_CUBE_QUERY_VALID_DIMENSIONS = [
    "gene_ontology_term_id",
    "tissue_ontology_term_id",
    "cell_type_ontology_term_id",
]
