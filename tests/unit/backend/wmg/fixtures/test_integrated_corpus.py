import pathlib
import time

import numpy as np
import tiledb
from scipy.sparse import csr_matrix, coo_matrix

from backend.wmg.data.schemas.corpus_schema import create_tdb, INTEGRATED_ARRAY_NAME


def create_test_corpus(corpus_name='test-corpus', snapshot_path=None, dataset_spec=([5, 3], [10, 5], [100, 30])):
    if not snapshot_path:
        timestamp = int(time.time())
        snapshot_path = f"{pathlib.Path().resolve()}/{timestamp}"
    corpus_path = f"{snapshot_path}/{corpus_name}"

    if not tiledb.VFS().is_dir(corpus_path):
        create_tdb(snapshot_path, corpus_name)

    for x in range(len(dataset_spec)):
        num_cells = x[0]
        num_genes= x[1]
        counts = csr_matrix(np.random.poisson(1, size=(num_cells, num_genes)), dtype=np.float32)
        update_corpus_obs()
        update_corpus_vars()
        update_corpus_counts()

    # consolidate corpus
    for arr_name in [f"{corpus_path}/{name}" for name in ["obs", "var", INTEGRATED_ARRAY_NAME]]:
        tiledb.consolidate(arr_name)
        tiledb.vacuum(arr_name)

    # assay_ontologies = np.random.choice(list(INCLUDED_ASSAYS.keys()), size=(adata.n_obs,))



def create_rankit_transformed integrated_matrix(corpus_path, gene_count=10):
    array_name = f"{corpus_path}/{INTEGRATED_ARRAY_NAME}"
    global_var_index = np.zeros((gene_count,), dtype=np.uint32)
    row = [0, 1, 2, 3]
    col = [0, 1, 2, 3]
    rankits = [0.3, 0.5, 0.7, 0.9]
    raw_counts = [4, 5, 8, 19]
    rankit_csr_matrix = csr_matrix((rankits, (row, col)))
    raw_counts_coo_matrix = coo_matrix((raw_counts, (row, col)))

    rankit_integrated_csr_matrix = rankit(raw_expression_csr_matrix)

rankit_integrated_coo_matrix = filter_out_rankits_with_low_expression_counts(
    rankit_integrated_csr_matrix, raw_expression_csr_matrix, expect_majority_filtered=True
)

global_rows = rankit_integrated_coo_matrix.row + start + first_obs_idx
global_cols = global_var_index[rankit_integrated_coo_matrix.col]

rankit_data = rankit_integrated_coo_matrix.data

    array[global_rows, global_cols] = {"rankit": rankit_data}
