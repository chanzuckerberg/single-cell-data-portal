import logging

import anndata
from scipy import sparse


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def validate_dataset_properties(anndata_object: anndata.AnnData) -> bool:
    if not sparse.issparse(anndata_object.X):
        logger.warning("No dense handling yet, not loading")
        return False
    schema_version = anndata_object.uns.get("schema_version", None)
    if not schema_version:
        logger.warning("Unknown schema, not loading")
        return False
    if schema_version < "2.0.0" or schema_version >= "3.0.0":
        logger.warning("Invalid schema version, not loading")
        return False
    return True


def validate_corpus_load(anndata_object: anndata.AnnData, group_name: str, dataset_id: str):
    """
    Validate that the load looks sane
    """
    pass
