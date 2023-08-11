import logging

import numba as nb
import numpy as np
import pandas as pd
import tiledb

from backend.common.utils.math_utils import MB
from backend.wmg.data.schemas.corpus_schema import INTEGRATED_ARRAY_NAME
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.data.utils import log_func_runtime
from backend.wmg.pipeline.summary_cubes.extract import extract_obs_data

logger = logging.getLogger(__name__)

BINOMIAL_NNZ_RANKIT_THR = 1


def transform(
    corpus_path: str, gene_ontology_term_ids: list, cube_dims: list
) -> (pd.DataFrame, np.ndarray, np.ndarray):  # type: ignore
    """
    Build the summary cube with rankit expression sum & sum of squares, nnz
    (num cells with non zero expression) values for each gene for each possible
    group of cell attributes (cube row). All values are aggregated across the
    descendants of the cell type in each group (row).
    """

    cell_labels, cube_index = make_cube_index(corpus_path, cube_dims)
    n_groups = len(cube_index)
    n_genes = len(gene_ontology_term_ids)

    cube_sum = np.zeros((n_groups, n_genes), dtype=np.float32)
    cube_sqsum = np.zeros((n_groups, n_genes), dtype=np.float32)
    cube_nnz = np.zeros((n_groups, n_genes), dtype=np.uint64)
    cube_nnz_thr = np.zeros((n_groups, n_genes), dtype=np.uint64)

    reduce_X(corpus_path, cell_labels.cube_idx.values, cube_sum, cube_sqsum, cube_nnz, cube_nnz_thr)
    return cube_index, cube_sum, cube_sqsum, cube_nnz, cube_nnz_thr


@log_func_runtime
def reduce_X(
    tdb_group: str,
    cube_indices: np.ndarray,
    cube_sum: np.ndarray,
    cube_sqsum: np.ndarray,
    cube_nnz: np.ndarray,
    cube_nnz_thr: np.ndarray,
):
    """
    Reduce the expression data stored in the integrated corpus by summing it by gene for each cube row (unique combo
    of cell attributes)
    """
    cfg = {
        "py.init_buffer_bytes": 512 * MB,
        "py.exact_init_buffer_bytes": "true",
    }
    with tiledb.open(f"{tdb_group}/{INTEGRATED_ARRAY_NAME}", ctx=create_ctx(config_overrides=cfg)) as expression:
        query_results = expression.query(return_incomplete=True, order="U", attrs=["rankit"])
        for i, result in enumerate(query_results.df[:]):
            logger.info(f"reduce integrated expression data, i={i}")
            gene_expression_sum_x_cube_dimension(
                result["rankit"].values,
                result["obs_idx"].values,
                result["var_idx"].values,
                cube_indices,
                cube_sum,
                cube_sqsum,
                cube_nnz,
                cube_nnz_thr,
            )


# TODO: this could be further optimize by parallel chunking.  Might help large arrays if compute ends up being a bottleneck. # noqa E501
@nb.njit(fastmath=True, error_model="numpy", parallel=False, nogil=True)
def gene_expression_sum_x_cube_dimension(
    rankit_values: np.ndarray,
    obs_idxs: np.ndarray,
    var_idx: np.ndarray,
    cube_indices: np.ndarray,
    sum_into: np.ndarray,
    sqsum_into: np.ndarray,
    nnz_into: np.ndarray,
    nnz_thr_into: np.ndarray,
):
    """
    Sum the rankit values for each gene (for each cube row/combo of cell attributes)
    Also track the number of cells that express that gene (nnz count)
    """
    for k in range(len(rankit_values)):
        val = rankit_values[k]
        if np.isfinite(val):
            cidx = var_idx[k]
            grp_idx = cube_indices[obs_idxs[k]]
            sum_into[grp_idx, cidx] += val
            sqsum_into[grp_idx, cidx] += val**2
            nnz_into[grp_idx, cidx] += 1
            nnz_thr_into[grp_idx, cidx] += val >= BINOMIAL_NNZ_RANKIT_THR


def make_cube_index(tdb_group: str, cube_dims: list) -> (pd.DataFrame, pd.DataFrame):  # type: ignore
    """
    Create index for queryable dimensions
    """
    cell_labels = extract_obs_data(tdb_group, cube_dims)

    # number of cells with specific tuple of dims
    cube_index = pd.DataFrame(cell_labels.value_counts(), columns=["n"])
    cube_index["cube_idx"] = range(len(cube_index))

    cell_labels = cell_labels.join(cube_index.cube_idx, on=cube_dims)
    # we failed to correctly create the corpus if these are false
    assert len(cell_labels.index) == cell_labels.index[-1] + 1
    assert cell_labels.index[0] == 0

    return cell_labels, cube_index
