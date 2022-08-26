import gc
import logging
import os
from typing import List

import tiledb

from backend.corpus_asset_pipelines.integrated_corpus import extract
from backend.corpus_asset_pipelines.integrated_corpus import load
from backend.corpus_asset_pipelines.integrated_corpus.transform import (
    apply_pre_concatenation_filters, create_high_level_tissue,
)
from backend.corpus_asset_pipelines.integrated_corpus.validate import should_load_dataset, validate_dataset_properties
from backend.wmg.data.schemas.corpus_schema import INTEGRATED_ARRAY_NAME, OBS_ARRAY_NAME, VAR_ARRAY_NAME
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.data.utils import log_func_runtime

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


@log_func_runtime
def extract_datasets(dataset_directory: List):
    s3_uris = extract.get_dataset_s3_uris()
    extract.copy_datasets_to_instance(s3_uris, dataset_directory)
    logger.info("Copied datasets to instance")


@log_func_runtime
def build_integrated_corpus(dataset_directory: List, corpus_path: str):
    """
    Given the path to a directory containing one or more h5ad files and a group name, call the h5ad loading function
    on all files, loading/concatenating the datasets together under the group name

    @mdunitz TODO: As the number of datasets grows this will become a bottleneck step.
                    hopefully this code will be removed in the near future as we switch to SOMA
                    if this becomes too slow before soma is ready for prod this step can be safely parallelized
                    with some consideration given to the load step to avoid collisions of the feature (obs) indices
                    in the integrated corpus
    """
    with tiledb.scope_ctx(create_ctx()):
        dataset_count = len(os.listdir(dataset_directory))
        for dataset in enumerate(os.listdir(dataset_directory)):
            logger.info(f"Processing dataset {dataset[0] + 1} of {dataset_count}")
            h5ad_file_path = f"{dataset_directory}/{dataset[1]}/local.h5ad"
            process_h5ad_for_corpus(
                h5ad_file_path, corpus_path
            )  # TODO Can this be parallelized? need to be careful handling global indexes but tiledb has a lock I think
            gc.collect()

        logger.info("all loaded, now consolidating.")

        for arr_name in [OBS_ARRAY_NAME, VAR_ARRAY_NAME, INTEGRATED_ARRAY_NAME]:
            arr_path = f"{corpus_path}/{arr_name}"
            tiledb.consolidate(arr_path)
            tiledb.vacuum(arr_path)


def process_h5ad_for_corpus(h5ad_path: str, corpus_path: str):
    """
    Given the location of a h5ad dataset and a group name, check the dataset is not already loaded
    then read the dataset into the tiledb object (under group name), updating the var and feature indexes
    to avoid collisions within the larger tiledb object
    """
    dataset_id = should_load_dataset(h5ad_path, corpus_path)
    if not dataset_id:
        return

    # extract
    anndata_object = extract.extract_h5ad(h5ad_path=h5ad_path)

    # transform
    apply_pre_concatenation_filters(anndata_object)
    create_high_level_tissue(anndata_object)
    logger.info(f"loaded: shape={anndata_object.shape}")
    if not validate_dataset_properties(anndata_object):
        return

    # load
    load.load_dataset(corpus_path, anndata_object, dataset_id)
