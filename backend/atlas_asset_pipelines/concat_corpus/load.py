import logging
from typing import List

import numpy as np
import pandas as pd
import tiledb

from backend.atlas_asset_pipelines.concat_corpus.transform import transform_dataset_raw_counts_to_rankit
from backend.wmg.data.schemas.corpus_schema import var_labels, obs_labels

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


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
