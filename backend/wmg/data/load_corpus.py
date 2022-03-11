import gc
import logging
import os
from typing import List, Union

import anndata
import numpy
import numpy as np
import pandas as pd
import tiledb
from anndata._core.views import ArrayView
from scipy import sparse

from backend.wmg.data.schemas.corpus_schema import var_labels, obs_labels

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def dataset_already_loaded(dataset_id: str, group_name: str) -> bool:
    if dataset_id in get_all_dataset_ids(group_name):
        logger.info("oops, that dataset is already loaded!")
        return True
    return False


def get_dataset_id(h5ad: str) -> str:
    dataset_id = os.path.splitext(os.path.split(h5ad)[1])[0]
    if dataset_id == 'local':
        dataset_id = os.path.split(os.path.split(h5ad)[0])[1]
    return dataset_id


def validate_dataset_properties(anndata_object: anndata.AnnData) -> bool:
    if not sparse.issparse(anndata_object.X):
        logger.warning("oops, no dense handling yet, not loading")
        return False
    if anndata_object.uns.get("schema_version", None) != "2.0.0":
        logger.warning("oops, unknown schema, not loading")
        return False
    return True


def load_h5ad(h5ad: str, group_name: str, validate: bool):
    """
    Given the location of a h5ad dataset and a group name, check the dataset is not already loaded
    then read the dataset into the tiledb object (under group name), updating the var and feature indexes
    to avoid collisions within the larger tiledb object
    """
    logger.info(f"Loading {h5ad}...")
    dataset_id = get_dataset_id(h5ad)
    if dataset_already_loaded(dataset_id, group_name):
        return

    anndata_object = anndata.read_h5ad(h5ad)
    logger.info(f"loaded: shape={anndata_object.shape}")
    if not validate_dataset_properties(anndata_object):
        return

    var_df = update_global_var(group_name, anndata_object.var)

    # Calculate mapping between var/feature coordinates in H5AD (file local) and TDB (global)
    global_var_index = np.zeros((anndata_object.shape[1],), dtype=np.uint32)
    var_feature_to_coord_map = {k: v for k, v in var_df[["gene_ontology_term_id", "var_idx"]].to_dict("split")["data"]}
    for idx in range(anndata_object.shape[1]):
        gene_ontology_term_id = anndata_object.var.index.values[idx]
        global_coord = var_feature_to_coord_map[gene_ontology_term_id]
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
        missing_var = set(src_var_df.index.to_numpy(dtype=str)) - set(
            var_df["gene_ontology_term_id"].to_numpy(dtype=str))

    if len(missing_var) > 0:
        logger.info(f"Adding {len(missing_var)} gene records...")
        missing_var_df = src_var_df[src_var_df.index.isin(missing_var)]
        save_axes_labels(missing_var_df, var_array_name, var_labels)
    with tiledb.open(var_array_name, "r") as var:
        var_df = var.df[:]
        var_df.index = var_df.gene_ontology_term_id

    logger.info(f"Global var index length: {var_df.shape}")
    return var_df


def save_axes_labels(df: pd.DataFrame, array_name: str, label_info: List) -> int:
    """
    # TODO
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


def save_X(anndata_object: anndata.AnnData, group_name: str, global_var_index: np.ndarray, first_obs_idx: int):
    """
    Save (pre)normalized expression counts to the tiledb corpus object
    """
    array_name = f"{group_name}/X"
    expression_matrix = anndata_object.X
    logger.debug(f"saving {array_name}...\n")
    stride = max(int(
        np.power(10, np.around(np.log10(1e9 / expression_matrix.shape[1])))
    ), 10_000)
    with tiledb.open(array_name, mode="w") as array:
        for start in range(0, expression_matrix.shape[0], stride):
            end = min(start + stride, expression_matrix.shape[0])
            sparse_expression_matrix = sparse.coo_matrix(expression_matrix[start:end, :])
            rows = sparse_expression_matrix.row + start + first_obs_idx
            cols = global_var_index[sparse_expression_matrix.col]
            data = sparse_expression_matrix.data

            array[rows, cols] = data
            del sparse_expression_matrix, rows, cols, data
            gc.collect()

    logger.debug(f"Saved {group_name}.")


def get_X_raw(anndata_object: anndata.AnnData) -> Union[np.ndarray, sparse.spmatrix, ArrayView]:
    """
    Current rules for our curated H5ADs:
    * if there is a .raw, it is the raw counts, and .X is transformed/normalized (by author) or is == to .raw.X
    * if there is no .raw, ad.X contains the raw counts
    """
    raw_expression_matrix = getattr(anndata_object.raw, "X", None)
    return raw_expression_matrix if raw_expression_matrix is not None else anndata_object.X


def save_raw(anndata_object: anndata.AnnData, group_name: str, global_var_index: numpy.ndarray, first_obs_idx: int):
    """
    Apply rankit normalization to raw expression values and save to the tiledb corpus object
    """
    array_name = f"{group_name}/raw"
    expression_matrix = get_X_raw(anndata_object)
    logger.info(f"saving {array_name}...")
    stride = max(int(
        np.power(10, np.around(np.log10(1e9 / expression_matrix.shape[1])))
    ), 10_000)
    with tiledb.open(array_name, mode="w") as array:
        for start in range(0, expression_matrix.shape[0], stride):
            end = min(start + stride, expression_matrix.shape[0])
            csr_sparse_raw_expression_matrix = sparse.csr_matrix(expression_matrix[start:end, :])

            coo_sparse_raw_expression_matrix = csr_sparse_raw_expression_matrix.tocoo(copy=False)
            rows = coo_sparse_raw_expression_matrix.row + start + first_obs_idx
            cols = global_var_index[coo_sparse_raw_expression_matrix.col]
            raw_data = coo_sparse_raw_expression_matrix.data

            rankit_normalized_coo_sparse_raw_expression_matrix = rankit(csr_sparse_raw_expression_matrix).tocoo(copy=False)
            assert np.array_equal(
                coo_sparse_raw_expression_matrix.row,
                rankit_normalized_coo_sparse_raw_expression_matrix.row
            )
            assert np.array_equal(
                coo_sparse_raw_expression_matrix.col,
                rankit_normalized_coo_sparse_raw_expression_matrix.col
            )
            rankit_data = rankit_normalized_coo_sparse_raw_expression_matrix.data

            array[rows, cols] = {"data": raw_data, "rankit": rankit_data}
            del coo_sparse_raw_expression_matrix, \
                rankit_normalized_coo_sparse_raw_expression_matrix, \
                rows, cols, raw_data, rankit_data
            gc.collect()

    logger.debug(f"Saved {array_name}.")
