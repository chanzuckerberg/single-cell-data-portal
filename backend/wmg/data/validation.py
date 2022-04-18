import gc
import logging

import numpy as np
import tiledb
from scipy import sparse

from backend.wmg.data.schemas.corpus_schema import obs_labels, create_local_to_global_gene_coord_index

logger = logging.getLogger(__name__)


def validate_corpus_load(ad, group_name, dataset_id):
    """
    Validate that the load matches base assumptions
    - no duplicates in gene table (var)
    - all obs_idx for a dataset are contiguous
    -
    """
    from backend.wmg.data.load_corpus import get_X_raw

    logger.info(f"validating...{group_name}, {dataset_id}")
    with tiledb.open(f"{group_name}/var") as var:
        with tiledb.open(f"{group_name}/obs") as obs:
            with tiledb.open(f"{group_name}/raw") as raw_X:

                logger.info("\tchecking var...")
                all_features = var.query().df[:]
                ## confirm no duplicates in gene table (var)
                assert all_features.shape[0] == len(set(all_features["gene_ontology_term_id"]))
                ## validate var
                assert all_features[all_features["gene_ontology_term_id"].isin(ad.var.index.values)].shape[0] == ad.n_vars

                ## validate obs
                logger.info("\tchecking obs...")
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
                logger.debug("\tchecking raw X...")
                var_idx_map = create_local_to_global_gene_coord_index(all_features, ad.var.index)
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


def validate_rankit_transformation():
    pass