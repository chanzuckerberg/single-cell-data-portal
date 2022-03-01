import logging
from typing import Union, List
import os.path
import gc

import numpy as np
import pandas as pd
from scipy import sparse
import tiledb

from .utils import get_all_dataset_ids, logmsg, create_fast_ctx
from .schema import obs_labels, var_labels
from .rankit import rankit

logger = logging.getLogger(__name__)


def load(dataset_files, group_name, validate):
    with tiledb.scope_ctx(create_fast_ctx()):
        if len(dataset_files) == 1 and os.path.isdir(dataset_files[0]):
            for dataset in os.listdir(dataset_files[0]):
                file_path = f"{dataset_files[0]}/{dataset}/local.h5ad"
                try:
                    load_h5ad(file_path, group_name, validate)
                except Exception as e:
                    logger.warning(f"Issue loading file: {dataset}, {e}")
                finally:
                    gc.collect()

        else:
            for dataset in dataset_files:
                try:
                    load_h5ad(dataset, group_name, validate)
                except Exception as e:
                    logger.warning(f"Issue loading file: {dataset}, {e}")
                finally:
                    gc.collect()

        logger.info("all loaded, now consolidating.")
        for arr_name in [f"{group_name}/{name}" for name in ["obs", "var", "raw", "X"]]:
            tiledb.consolidate(arr_name)
            tiledb.vacuum(arr_name)


def load_h5ad(h5ad, group_name, validate):
    import anndata

    logger.info(f"Loading {h5ad}...")

    dataset_id = os.path.splitext(os.path.split(h5ad)[1])[0]
    if dataset_id == 'local':
        dataset_id = os.path.split(os.path.split(h5ad)[0])[1]
    if dataset_id in get_all_dataset_ids(group_name):
        logger.info("oops, that dataset is already loaded!")
        return

    ad = anndata.read_h5ad(h5ad)
    logger.info(f"loaded: shape={ad.shape}")
    if not sparse.issparse(ad.X):
        logger.warning("oops, no dense handling yet, not loading")
        return
    if ad.uns.get("schema_version", None) != "2.0.0":
        logger.warning("oops, unknown schema, not loading")
        return

    var_df = update_global_var(group_name, ad.var)

    # Calculate mapping between var/feature coordinates in H5AD (file local) and TDB (global)
    global_var_index = np.zeros((ad.shape[1],), dtype=np.uint32)
    var_feature_to_coord_map = {k: v for k, v in var_df[["feature_id", "var_idx"]].to_dict("split")["data"]}
    for idx in range(ad.shape[1]):
        feature_id = ad.var.index.values[idx]
        global_coord = var_feature_to_coord_map[feature_id]
        global_var_index[idx] = global_coord

    obs = ad.obs
    obs["dataset_id"] = dataset_id
    first_obs_idx = save_axes_labels(obs, f"{group_name}/obs", obs_labels)
    save_raw(ad, group_name, global_var_index, first_obs_idx)
    save_X(ad, group_name, global_var_index, first_obs_idx)

    if validate:
        validate_load(ad, group_name, dataset_id)


def update_global_var(group_name: str, src_var_df: pd.DataFrame) -> pd.DataFrame:
    """
    Update the global var (feature) array. Adds any missing features we have
    not seen before.

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

    logger.info("\tglobal var index length: ", var_df.shape)
    return var_df


def create_local_to_global_feature_coord_index(
        var_df: pd.DataFrame, feature_ids: Union[List[str], np.ndarray]
) -> np.ndarray:
    """
    Create an array mapping feature ids local to global index
    """
    n_features = len(feature_ids)
    local_to_global_feature_coord = np.zeros((n_features,), dtype=np.uint32)
    var_feature_to_coord_map = {k: v for k, v in var_df[["feature_id", "var_idx"]].to_dict("split")["data"]}
    for idx in range(n_features):
        feature_id = feature_ids[idx]
        global_coord = var_feature_to_coord_map[feature_id]
        local_to_global_feature_coord[idx] = global_coord

    return local_to_global_feature_coord


def save_axes_labels(df, array_name, label_info) -> int:
    logger.info(f"saving {array_name}...", end="", flush=True)

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
    array_name = f"{group_name}/raw"
    X = get_X_raw(ad)
    logger.info(f"saving {array_name}...", end="", flush=True)
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
    array_name = f"{group_name}/X"
    X = ad.X
    logger.info(f"saving {array_name}...", end="", flush=True)
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


def validate_load(ad, group_name, dataset_id):
    """
    Validate that the load looks sane
    """

    logmsg("validating...")
    with tiledb.open(f"{group_name}/var") as var:
        with tiledb.open(f"{group_name}/obs") as obs:
            with tiledb.open(f"{group_name}/raw") as raw_X:

                print("\tchecking var...")
                all_features = var.query().df[:]
                ## confirm no duplicates in gene table (var)
                assert all_features.shape[0] == len(set(all_features["feature_id"]))
                ## validate var
                assert all_features[all_features["feature_id"].isin(ad.var.index.values)].shape[0] == ad.n_vars

                ## validate obs
                print("\tchecking obs...")
                tdb_obs = obs.df[dataset_id].sort_values(by=["obs_idx"], ignore_index=True)
                start_coord = tdb_obs.iloc[0].obs_idx
                assert start_coord == tdb_obs.obs_idx.min()
                h5ad_df = ad.obs
                h5ad_df["dataset_id"] = dataset_id
                assert tdb_obs.shape[0] == h5ad_df.shape[0]
                for lbl in obs_labels:
                    datum = lbl.decode(h5ad_df, start_coord=start_coord)
                    tdb_data = tdb_obs[lbl.key].values
                    if lbl.dtype in ["ascii", np.bytes_]:
                        # tiledb .df converts np.bytes_ back to str in some cases.  Arg.
                        datum = datum.astype(str)
                        tdb_data = tdb_data.astype(str)
                    assert (datum == tdb_data).all()

                # We assume that all obs_idx for a dataset are contiguous.  Useful assumption.
                obs_idx = obs.df[dataset_id].obs_idx
                assert obs_idx.max() - obs_idx.min() + 1 == obs_idx.shape[0]

                ## Validate X
                print("\tchecking raw X...")
                var_idx_map = create_local_to_global_feature_coord_index(all_features, ad.var.index)
                stride = 100_000
                starting_obs_idx = obs.df[dataset_id].obs_idx.min()
                for start in range(0, ad.n_obs, stride):
                    end = min(start + stride, ad.n_obs)

                    # from H5AD - must map column indices into the global space to compare
                    adX = get_X_raw(ad)
                    h5ad_X_slice_coo = adX[start:end, :].tocoo()
                    h5ad_X_slice_remapped = sparse.coo_matrix(
                        (h5ad_X_slice_coo.data, (h5ad_X_slice_coo.row, var_idx_map[h5ad_X_slice_coo.col]))
                    ).tocsr()
                    del h5ad_X_slice_coo

                    # from TDB
                    tdb_X_slice = raw_X.query(attrs=["data"])[(starting_obs_idx + start): (starting_obs_idx + end), :]
                    tdb_X_slice_remapped = sparse.coo_matrix(
                        (
                            tdb_X_slice["data"],
                            (tdb_X_slice["obs_idx"] - start - starting_obs_idx, tdb_X_slice["var_idx"]),
                        )
                    ).tocsr()
                    del tdb_X_slice

                    assert (h5ad_X_slice_remapped.shape == tdb_X_slice_remapped.shape) and (
                            h5ad_X_slice_remapped != tdb_X_slice_remapped
                    ).nnz == 0

                    del h5ad_X_slice_remapped, tdb_X_slice_remapped
                    gc.collect()


def roundHalfToEven(a: np.ndarray, keepbits: int) -> np.ndarray:
    """
    Generate reduced precision floating point array.

    Ref: https://gmd.copernicus.org/articles/14/377/2021/gmd-14-377-2021.html
    """
    assert a.dtype == np.float32
    if keepbits < 1 or keepbits >= 23:
        return a
    b = a.view(dtype=np.int32)
    maskbits = 23 - keepbits
    mask = (0xFFFFFFFF >> maskbits) << maskbits
    half_quantum1 = (1 << (maskbits - 1)) - 1
    b += ((b >> maskbits) & 1) + half_quantum1
    b &= mask
    return a
