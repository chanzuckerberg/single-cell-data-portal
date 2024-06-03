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
