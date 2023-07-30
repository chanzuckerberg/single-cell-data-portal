# When the API is set to read a particular
# snapshot schema version, it will load the
# latest snapshot for the schema version by
# default.
#
# In the case a specfic snapshot must be loaded
# WMG_API_FORCE_LOAD_SNAPSHOT_ID should be set.
WMG_API_SNAPSHOT_SCHEMA_VERSION = "v1"

# In the case we need to rollback or rollforward
# set this variable to a specific snapshot id
# that should be loaded.
#
# NOTE: The snapshot id that should be force
# loaded must belong to the schema version set
# in WMG_API_SNAPSHOT_SCHEMA_VERSION
WMG_API_FORCE_LOAD_SNAPSHOT_ID = None
