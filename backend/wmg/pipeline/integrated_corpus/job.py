import gc
import json
import logging
import os
from typing import List, Tuple, Union

import tiledb

from backend.wmg.data.schemas.corpus_schema import (
    DATASET_TO_GENE_IDS_NAME,
    INTEGRATED_ARRAY_NAME,
    OBS_ARRAY_NAME,
    VAR_ARRAY_NAME,
)
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.data.utils import log_func_runtime
from backend.wmg.pipeline.integrated_corpus import extract, load
from backend.wmg.pipeline.integrated_corpus.transform import apply_pre_concatenation_filters, create_high_level_tissue
from backend.wmg.pipeline.integrated_corpus.validate import should_load_dataset, validate_dataset_properties
from backend.wmg.pipeline.utils import remap_anndata_normalized_X_to_raw_X_if_exists

logger: logging.Logger = logging.getLogger(__name__)


@log_func_runtime
def extract_datasets(dataset_directory: List):
    asset_urls = extract.get_dataset_asset_urls()
    extract.copy_datasets_to_instance(asset_urls, dataset_directory)
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
        if os.path.exists(f"{corpus_path}/{DATASET_TO_GENE_IDS_NAME}.json"):
            with open(f"{corpus_path}/{DATASET_TO_GENE_IDS_NAME}.json", "r") as fp:
                dataset_gene_mapping = json.load(fp)
        else:
            dataset_gene_mapping = {}

        for index, dataset in enumerate(os.listdir(dataset_directory)):
            logger.info(f"Processing dataset {index + 1} of {dataset_count}")
            h5ad_file_path = f"{dataset_directory}/{dataset}/local.h5ad"
            logger.info(f"{h5ad_file_path=}")
            dataset_id, gene_ids = process_h5ad_for_corpus(h5ad_file_path, corpus_path)
            # (mdunitz) TODO Can the above be parallelized? need to be careful handling
            # global indexes but tiledb has a lock I think
            gc.collect()
            if dataset_id and gene_ids:
                dataset_gene_mapping[dataset_id] = gene_ids

        with tiledb.open(f"{corpus_path}/{VAR_ARRAY_NAME}") as var:
            gene_count = len(var.query().df[:])
        with tiledb.open(f"{corpus_path}/{OBS_ARRAY_NAME}") as obs:
            cell_count = len(obs.query().df[:])
        with open(f"{corpus_path}/{DATASET_TO_GENE_IDS_NAME}.json", "w") as d2g:
            json.dump(dataset_gene_mapping, d2g)
        logger.info("all loaded, now consolidating.")
        for arr_name in [OBS_ARRAY_NAME, VAR_ARRAY_NAME, INTEGRATED_ARRAY_NAME]:
            arr_path = f"{corpus_path}/{arr_name}"
            tiledb.consolidate(arr_path)
            tiledb.vacuum(arr_path)

    logger.info(f"{dataset_count=}, {gene_count=}, {cell_count=}")


def process_h5ad_for_corpus(h5ad_path: str, corpus_path: str) -> Union[Tuple[str, dict], Tuple[None, None]]:
    """
    Given the location of a h5ad dataset and a group name, check the dataset is not already loaded
    then read the dataset into the tiledb object (under group name), updating the var and feature indexes
    to avoid collisions within the larger tiledb object

    Returns a tuple of (dataset_id: str, gene_ids: list) if the dataset is loaded and validated. Otherwise,
    returns a tuple of (None, None).
    """
    dataset_id = should_load_dataset(h5ad_path, corpus_path)
    if not dataset_id:
        return None, None

    # extract
    try:
        anndata_object = extract.extract_h5ad(h5ad_path=h5ad_path)
    except Exception as e:
        logger.error(f"Failed to extract h5ad file {h5ad_path} with error: {e}")
        return None, None

    # A memory optimization to reduce memory consumption when loading large datasets
    # is to remap AnnData.X = AnnData.raw.X if AnnData.raw.X exists and then garbage
    # collect the original object pointed to by AnnData.X.
    #
    # This is significant memory savings because size_of(AnnData.X) ~= size_of(AnnData.raw.X)
    # and so we reduce memory pressure by 50%. Moreover, the normalized value in the
    # original AnnData.X is not used throughout the pipeline run.
    remap_anndata_normalized_X_to_raw_X_if_exists(anndata_object)

    # transform
    apply_pre_concatenation_filters(anndata_object)
    create_high_level_tissue(anndata_object)
    logger.info(f"loaded: shape={anndata_object.shape}")
    if not validate_dataset_properties(anndata_object):
        return None, None

    # load
    load.load_dataset(corpus_path, anndata_object, dataset_id)
    return dataset_id, list(anndata_object.var_names)
