# When this config flag is set, the API will load the snapshot
# from the local disk. When the flag is False, the API
# will load the snapshot from S3
CENSUS_CUBE_API_READ_FS_CACHED_SNAPSHOT = True

# When the API is set to read a particular
# snapshot schema version, it will load the
# latest snapshot for the schema version by
# default.
#
# In the case a specfic snapshot must be loaded
# CENSUS_CUBE_API_FORCE_LOAD_SNAPSHOT_ID should be set.
# LATEST_READER_SNAPSHOT_SCHEMA_VERSION in container_init.sh and CENSUS_CUBE_API_SNAPSHOT_SCHEMA_VERSION
# below should have the same value.
CENSUS_CUBE_API_SNAPSHOT_SCHEMA_VERSION = "v5"

# In the case we need to rollback or rollforward
# set this variable to a specific snapshot id
# that should be loaded.
#
# NOTE: The snapshot id that should be force
# loaded must belong to the schema version set
# in CENSUS_CUBE_API_SNAPSHOT_SCHEMA_VERSION
CENSUS_CUBE_API_FORCE_LOAD_SNAPSHOT_ID = None


# These are the valid attributes and dimensions consulted by the
# wmg api (reader) to determine the list of attributes and dimensions
# TO RETRIEVE from the cube when performing a query.
# This is important because the cube allocates a large amount of memory for
# each attribute and dimension retrieved. Therefore, strictly specifying the
# attributes and dimensions to retrieve makes the query more memory efficient
READER_CENSUS_CUBE_CUBE_QUERY_VALID_ATTRIBUTES = ["gene_ontology_term_id", "cell_type_ontology_term_id"]
READER_CENSUS_CUBE_CUBE_QUERY_VALID_DIMENSIONS = [
    "gene_ontology_term_id",
    "tissue_ontology_term_id",
]
