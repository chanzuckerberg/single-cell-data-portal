class CorporaConstants:
    REQUIRED_SIMPLE_METADATA_FIELDS = [
        "schema_version",
        "title",
    ]

    # Lists are encoded as ndarrays in h5ad/anndata, so must be cast to lists during cxg conversion
    OPTIONAL_LIST_METADATA_FIELDS = [
        "batch_condition",
    ]

    OPTIONAL_SIMPLE_METADATA_FIELDS = [
        "default_embedding",
        "X_approximate_distribution",
    ]
    SUPER_CURATOR_NAME = "super-curator"
    SUPER_CURATOR_SCOPE = "write:collections"
    SUPER_CURATOR_SUBMISSON_PATH = "super"

    CXG_ADMIN_SCOPE = "delete:collections"

    ORIGINAL_H5AD_ARTIFACT_FILENAME = "raw.h5ad"
    LABELED_H5AD_ARTIFACT_FILENAME = "local.h5ad"
