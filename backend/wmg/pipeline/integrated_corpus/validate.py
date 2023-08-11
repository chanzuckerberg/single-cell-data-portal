import logging
import os
from functools import cache

import anndata
from scipy import sparse

from backend.wmg.data.utils import get_all_dataset_ids
from backend.wmg.pipeline.integrated_corpus.extract import get_X_raw

logger = logging.getLogger(__name__)


def should_load_dataset(h5ad_path: str, corpus_path: str) -> str:
    dataset_id = get_dataset_id(h5ad_path)
    if is_dataset_already_loaded(corpus_path, dataset_id):
        return None  # type: ignore
    return dataset_id


@cache
def is_dataset_already_loaded(corpus_path: str, dataset_id: str) -> bool:
    if dataset_id in get_all_dataset_ids(corpus_path):
        logger.info(f"oops, {dataset_id=} is already loaded!")
        return True
    return False


def get_dataset_id(h5ad_path: str) -> str:
    dataset_id = os.path.splitext(os.path.split(h5ad_path)[1])[0]
    if dataset_id == "local":
        dataset_id = os.path.split(os.path.split(h5ad_path)[0])[1]
    return dataset_id


def validate_dataset_properties(anndata_object: anndata.AnnData) -> bool:
    expression_matrix = get_X_raw(anndata_object)
    if not sparse.issparse(expression_matrix):
        logger.warning("No dense handling yet, not loading")
        return False
    schema_version = anndata_object.uns.get("schema_version", None)
    if not schema_version:
        logger.warning("Unknown schema, not loading")
        return False
    if schema_version < "3.0.0":
        logger.warning("Invalid schema version, not loading")
        return False
    return True
