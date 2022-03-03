import logging
from typing import List
import os.path
import gc

import numpy as np
import pandas as pd
from scipy import sparse
import tiledb

from .utils import get_all_dataset_ids
from backend.wmg.data.config import create_fast_ctx
from .corpus_schema import obs_labels, var_labels
from .rankit import rankit
from .validate import validate_load

logger = logging.getLogger(__name__)


def load(dataset_directory: List, group_name: str, validate: bool):
    """
    Given the path to a directory containing one or more h5ad files and a group name, call the cube loading function
    on all files, loading/concatenating the datasets together under the group name
    """
    with tiledb.scope_ctx(create_fast_ctx()):
        for dataset in os.listdir(dataset_directory):
            file_path = f"{dataset_directory}/{dataset}/local.h5ad"
            try:
                load_h5ad(file_path, group_name, validate) # TODO Can this be parallelized? need to be careful handling global indexes
            except Exception as e:
                logger.warning(f"Issue loading file: {dataset}, {e}")
            finally:
                gc.collect()

        logger.info("all loaded, now consolidating.")
        for arr_name in [f"{group_name}/{name}" for name in ["obs", "var", "raw", "X"]]:
            tiledb.consolidate(arr_name)
            tiledb.vacuum(arr_name)


def load_h5ad(h5ad, group_name, validate):
    """
    Given the location of a h5ad dataset and a group name, check the dataset is not already loaded 
    then read the dataset into the tiledb object (under group name), updating the var and feature indexes 
    to avoid collisions within the larger tiledb object
    """
    import anndata
    logger.info(f"Loading {h5ad}...")

    dataset_id = os.path.splitext(os.path.split(h5ad)[1])[0]
    if dataset_id == 'local':
        dataset_id = os.path.split(os.path.split(h5ad)[0])[1]
    if dataset_id in get_all_dataset_ids(group_name):
        logger.info("oops, that dataset is already loaded!")
        return

    anndata_object = anndata.read_h5ad(h5ad)
    logger.info(f"loaded: shape={anndata_object.shape}")
    if not sparse.issparse(anndata_object.X):
        logger.warning("oops, no dense handling yet, not loading")
        return
    if anndata_object.uns.get("schema_version", None) != "2.0.0":
        logger.warning("oops, unknown schema, not loading")
        return

    var_df = update_global_var(group_name, anndata_object.var)

    # Calculate mapping between var/feature coordinates in H5AD (file local) and TDB (global)
    global_var_index = np.zeros((anndata_object.shape[1],), dtype=np.uint32)
    var_feature_to_coord_map = {k: v for k, v in var_df[["feature_id", "var_idx"]].to_dict("split")["data"]}
    for idx in range(anndata_object.shape[1]):
        feature_id = anndata_object.var.index.values[idx]
        global_coord = var_feature_to_coord_map[feature_id]
        global_var_index[idx] = global_coord

    obs = anndata_object.obs
    obs["dataset_id"] = dataset_id
    first_obs_idx = save_axes_labels(obs, f"{group_name}/obs", obs_labels)
    save_raw(anndata_object, group_name, global_var_index, first_obs_idx)
    save_X(anndata_object, group_name, global_var_index, first_obs_idx)

    if validate:
        validate_load(anndata_object, group_name, dataset_id)


def update_global_var(group_name: str, src_var_df: pd.DataFrame) -> pd.DataFrame:
    """
    Update the global var (gene) array. Adds any gene_ids we have not seen before.
    Returns the global var array as dataframe
    """
    var_array_name = f"{group_name}/var"
    with tiledb.open(var_array_name, "r") as var:
        var_df = var.df[:]
        missing_var = set(src_var_df.index.to_numpy(dtype=str)) - set(var_df["feature_id"].to_numpy(dtype=str))

    if len(missing_var) > 0:
        logger.info(f"Adding {len(missing_var)} gene records...")
        missing_var_df = src_var_df[src_var_df.index.isin(missing_var)]
        save_axes_labels(missing_var_df, var_array_name, var_labels)

    with tiledb.open(var_array_name, "r") as var:
        var_df = var.df[:]
        var_df.index = var_df.feature_id

    logger.info(f"Global var index length: {var_df.shape}")
    return var_df


def save_axes_labels(df, array_name, label_info) -> int:
    """
    # TODO
    """
    logger.info(f"saving {array_name}...")

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


def get_X_raw(ad):
    """
    Current rules for our curated H5ADs:
    * if there is a .raw, it is the raw counts, and .X is transformed/normalized (by author) or is == to .raw.X
    * if there is no .raw, ad.X contains the raw counts
    """
    raw_X = getattr(ad.raw, "X", None)
    return raw_X if raw_X is not None else ad.X


def save_raw(ad, group_name, global_var_index, first_obs_idx):
    """
    Apply rankit normalization to raw expression values and save to the tiledb corpus object
    """
    array_name = f"{group_name}/raw"
    X = get_X_raw(ad)
    logger.info(f"saving {array_name}...")
    stride = max(int(np.power(10, np.around(np.log10(1e9 / X.shape[1])))), 10_000)
    with tiledb.open(array_name, mode="w") as array:
        for start in range(0, X.shape[0], stride):
            end = min(start + stride, X.shape[0])
            X_raw = sparse.csr_matrix(X[start:end, :])
            X_raw_coo = X_raw.tocoo(copy=False)
            rows = X_raw_coo.row + start + first_obs_idx
            cols = global_var_index[X_raw_coo.col]
            raw_data = X_raw_coo.data

            X_rankit = rankit(X_raw)
            X_rankit_coo = X_rankit.tocoo(copy=False)
            assert np.array_equal(X_raw_coo.row, X_rankit_coo.row)
            assert np.array_equal(X_raw_coo.col, X_rankit_coo.col)
            rankit_data = X_rankit_coo.data

            array[rows, cols] = {"data": raw_data, "rankit": rankit_data}
            del X_raw_coo, X_rankit_coo, rows, cols, raw_data, rankit_data
            gc.collect()

    logger.info("saved.")


def save_X(ad, group_name, global_var_index, first_obs_idx):
    """
    Save (pre)normalized expression counts to the tiledb corpus object
    """
    array_name = f"{group_name}/X"
    X = ad.X
    logger.info(f"saving {array_name}...")
    stride = max(int(np.power(10, np.around(np.log10(1e9 / X.shape[1])))), 10_000)
    with tiledb.open(array_name, mode="w") as array:
        for start in range(0, X.shape[0], stride):
            end = min(start + stride, X.shape[0])
            X_coo = sparse.coo_matrix(X[start:end, :])
            rows = X_coo.row + start + first_obs_idx
            cols = global_var_index[X_coo.col]
            data = X_coo.data

            array[rows, cols] = data
            del X_coo, rows, cols, data
            gc.collect()

    logger.info("saved.")

