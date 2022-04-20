import concurrent
import logging
import time

import numba as nb
import numpy as np
import pandas as pd

import tiledb

from backend.atlas_asset_pipelines.expression_summary_cube.extract import extract_obs_data
from backend.corpora.common.utils.math_utils import MB
from backend.wmg.data.schemas.corpus_schema import INTEGRATED_ARRAY_NAME
from backend.wmg.data.schemas.cube_schema import cube_non_indexed_dims
from backend.wmg.data.tiledb import create_ctx

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


cube_indexed_dims_no_gene_ontology = [
    "tissue_ontology_term_id",
    "organism_ontology_term_id",
]


# todo name better
def transform(corpus_path: str, gene_ontology_term_ids: list):
    start_time = time.time()

    ##
    # Reduce X
    ##
    cube_dims = cube_indexed_dims_no_gene_ontology + cube_non_indexed_dims
    cell_labels, cube_index = make_cube_index(corpus_path, cube_dims)
    n_groups = len(cube_index)
    n_genes = len(gene_ontology_term_ids)

    cube_sum = np.zeros((n_groups, n_genes), dtype=np.float32)
    cube_nnz = np.zeros((n_groups, n_genes), dtype=np.uint64)
    cube_min = np.zeros((n_groups, n_genes), dtype=np.float32)
    cube_max = np.zeros((n_groups, n_genes), dtype=np.float32)

    # pass 1 - sum, nnz, min, max
    reduce_X(corpus_path, start_time, cell_labels.cube_idx.values, cube_sum, cube_nnz, cube_min, cube_max)
    return cube_index, cube_sum, cube_nnz


def reduce_X(tdb_group, start_time, cube_indices, *accum):
    """
    # TODO
    """
    with concurrent.futures.ThreadPoolExecutor() as tp:
        cfg = {
            "py.init_buffer_bytes": 512 * MB,
            "py.exact_init_buffer_bytes": "true",
        }
        with tiledb.open(f"{tdb_group}/{INTEGRATED_ARRAY_NAME}", ctx=create_ctx(config_overrides=cfg)) as X:
            iterable = X.query(return_incomplete=True, order="U", attrs=["rankit"])
            future = None
            for i, result in enumerate(iterable.df[:]):
                logger.info(f"reduce integrated expression data, iter {i}, {time.time() - start_time}")
                if future is not None:
                    future.result()  # forces a wait
                future = tp.submit(
                    coo_cube_pass1_into,
                    result["rankit"].values,
                    result["obs_idx"].values,
                    result["var_idx"].values,
                    cube_indices,
                    *accum,
                )

        return accum


"""
TODO: this could be further optimize by parallel chunking.  Might help large
arrays if compute ends up being a bottleneck.
"""


@nb.njit(fastmath=True, error_model="numpy", parallel=False, nogil=True)
def coo_cube_pass1_into(data, row, col, row_groups, sum_into, nnz_into, min_into, max_into):
    """
    # TODO
    """
    for k in range(len(data)):
        val = data[k]
        if np.isfinite(val):
            cidx = col[k]
            grp_idx = row_groups[row[k]]
            sum_into[grp_idx, cidx] += val
            nnz_into[grp_idx, cidx] += 1
            if val < min_into[grp_idx, cidx]:
                min_into[grp_idx, cidx] = val
            if val > max_into[grp_idx, cidx]:
                max_into[grp_idx, cidx] = val


def make_cube_index(tdb_group, cube_dims):
    """
    Create index for queryable dimensions
    """
    cell_labels = extract_obs_data(tdb_group, cube_dims)

    # number of cells with specific tuple of dims
    cube_index = pd.DataFrame(cell_labels.value_counts(), columns=["n"])
    cube_index["cube_idx"] = range(len(cube_index))

    #
    cell_labels = cell_labels.join(cube_index.cube_idx, on=cube_dims)

    # we failed to correctly create the cube if these are false
    assert len(cell_labels.index) == cell_labels.index[-1] + 1
    assert cell_labels.index[0] == 0

    return cell_labels, cube_index
