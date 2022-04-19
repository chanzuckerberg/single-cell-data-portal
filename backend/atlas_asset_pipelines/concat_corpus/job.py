import gc
import logging
import os
from typing import List

import anndata
import tiledb

from backend.atlas_asset_pipelines.concat_corpus.load import load_h5ad
from backend.atlas_asset_pipelines.concat_corpus.transform import apply_pre_concatenation_filters
from backend.atlas_asset_pipelines.concat_corpus.validate import validate_dataset_properties, validate_corpus_load
from backend.wmg.data.schemas.corpus_schema import INTEGRATED_ARRAY_NAME
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.data.utils import get_all_dataset_ids

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def build_concat_corpus(dataset_directory: List, corpus_path: str, validate: bool = False):
    """
    Given the path to a directory containing one or more h5ad files and a group name, call the h5ad loading function
    on all files, loading/concatenating the datasets together under the group name
    """
    with tiledb.scope_ctx(create_ctx()):
        dataset_count = len(os.listdir(dataset_directory))
        i = 0
        for dataset in os.listdir(dataset_directory):
            i += 1
            logger.info(f"Processing dataset {i} of {dataset_count}")
            h5ad_file_path = f"{dataset_directory}/{dataset}/local.h5ad"
            process_h5ad_for_corpus(
                h5ad_file_path, corpus_path, validate
            )  # TODO Can this be parallelized? need to be careful handling global indexes but tiledb has a lock I think
            gc.collect()

        logger.info("all loaded, now consolidating.")
        for arr_name in [f"{corpus_path}/{name}" for name in ["obs", "var", INTEGRATED_ARRAY_NAME]]:
            tiledb.consolidate(arr_name)
            tiledb.vacuum(arr_name)


def process_h5ad_for_corpus(h5ad_path: str, corpus_path: str, validate: bool):
    """
    Given the location of a h5ad dataset and a group name, check the dataset is not already loaded
    then read the dataset into the tiledb object (under group name), updating the var and feature indexes
    to avoid collisions within the larger tiledb object
    """
    logger.info(f"Loading {h5ad_path}...")
    dataset_id = get_dataset_id(h5ad_path)
    if is_dataset_already_loaded(corpus_path, dataset_id):
        return

    # extract
    anndata_object = anndata.read_h5ad(h5ad_path)

    # transform
    apply_pre_concatenation_filters(anndata_object)

    logger.info(f"loaded: shape={anndata_object.shape}")
    if not validate_dataset_properties(anndata_object):
        return

    # load
    load_h5ad(corpus_path, anndata_object, dataset_id)

    if validate:
        validate_corpus_load(anndata_object, corpus_path, dataset_id)


def is_dataset_already_loaded(corpus_path: str, dataset_id: str) -> bool:
    if dataset_id in get_all_dataset_ids(corpus_path):
        logger.info("oops, that dataset is already loaded!")
        return True
    return False


def get_dataset_id(h5ad_path: str) -> str:
    dataset_id = os.path.splitext(os.path.split(h5ad_path)[1])[0]
    if dataset_id == "local":
        dataset_id = os.path.split(os.path.split(h5ad_path)[0])[1]
    return dataset_id
