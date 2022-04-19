import gc
import logging
import os
from typing import List

import anndata
import numpy as np
import pandas as pd
import tiledb

from backend.atlas_asset_pipelines.concat_corpus.transform import transform_dataset_raw_counts_to_rankit, \
    apply_pre_concatenation_filters
from backend.atlas_asset_pipelines.concat_corpus.validate import validate_dataset_properties
from backend.wmg.data.schemas.corpus_schema import INTEGRATED_ARRAY_NAME, var_labels, obs_labels
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.data.utils import get_all_dataset_ids
from backend.wmg.data.validation import validate_corpus_load

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def load_h5ad_datasets(dataset_directory: List, corpus_path: str, validate: bool = False):
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


def load_h5ad(corpus_path, anndata_object, dataset_id):
    var_df = update_global_var(corpus_path, anndata_object.var)

    # Calculate mapping between var/feature coordinates in H5AD (file local) and TDB (global)
    # todo use gene_ontology_id as var_idx
    global_var_index = np.zeros((anndata_object.shape[1],), dtype=np.uint32)
    var_feature_to_coord_map = {k: v for k, v in var_df[["gene_ontology_term_id", "var_idx"]].to_dict("split")["data"]}
    for idx in range(anndata_object.shape[1]):
        gene_ontology_term_id = anndata_object.var.index.values[idx]
        global_coord = var_feature_to_coord_map[gene_ontology_term_id]
        global_var_index[idx] = global_coord

    obs = anndata_object.obs
    obs["dataset_id"] = dataset_id
    first_obs_idx = save_axes_labels(obs, f"{corpus_path}/obs", obs_labels)
    transform_dataset_raw_counts_to_rankit(anndata_object, corpus_path, global_var_index, first_obs_idx)


# todo rename
def save_axes_labels(df: pd.DataFrame, array_name: str, label_info: List) -> int:
    """
    Safely add given dataframe to extend the specified array with the appropriate encoding
    typically used to update concat_corpus obs/var frames for each new dataset

    """
    logger.info(f"Saving {array_name}...\n")

    with tiledb.open(array_name) as array:
        next_join_index = array.meta.get("next_join_index", 0)

    with tiledb.open(array_name, mode="w") as array:
        data = {}
        coords = []
        for lbl in label_info:
            datum = lbl.decode(df, next_join_index)
            if lbl.encode_as_dim:
                coords.append(datum)
            else:
                data[lbl.key] = datum
        array[tuple(coords)] = data
        array.meta["next_join_index"] = next_join_index + len(coords[0])

    logger.info("saved.")
    return next_join_index


def update_global_var(corpus_path: str, src_var_df: pd.DataFrame) -> pd.DataFrame:
    """
    Update the global var (gene) array. Adds any gene_ids we have not seen before.
    Returns the global var array as dataframe
    """

    var_array_name = f"{corpus_path}/var"
    with tiledb.open(var_array_name, "r") as var:
        var_df = var.df[:]
        missing_var = set(src_var_df.index.to_numpy(dtype=str)) - set(
            var_df["gene_ontology_term_id"].to_numpy(dtype=str)
        )

    if len(missing_var) > 0:
        logger.info(f"Adding {len(missing_var)} gene records...")
        missing_var_df = src_var_df[src_var_df.index.isin(missing_var)]
        save_axes_labels(missing_var_df, var_array_name, var_labels)
    with tiledb.open(var_array_name, "r") as var:
        var_df = var.df[:]
        var_df.index = var_df.gene_ontology_term_id

    logger.info(f"Global var index length: {var_df.shape}")
    return var_df


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

# TODO where does this belong
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