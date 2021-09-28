

class CorporaConstants(object):
    REQUIRED_SIMPLE_METADATA_FIELDS = [
        "schema_version",
        "title",
        "X_normalization",
    ]

    # Lists are encoded as ndarrays in h5ad/anndata, so must be cast to lists during cxg conversion
    OPTIONAL_LIST_METADATA_FIELDS = [
        "batch_condition",
    ]

    OPTIONAL_SIMPLE_METADATA_FIELDS = [
        "default_embedding",
        "X_approximate_distribution",
    ]
