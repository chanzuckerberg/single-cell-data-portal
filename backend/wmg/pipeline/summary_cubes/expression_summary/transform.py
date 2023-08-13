import logging
from dataclasses import dataclass
from typing import Tuple

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


@dataclass
class TransformResult:
    cube_index: pd.DataFrame
    cube_sum: np.ndarray
    cube_nnz: np.ndarray
    cube_sqsum: np.ndarray


def transform(*, corpus_path: str, gene_ontology_term_ids: list, cube_dims: list) -> TransformResult:
    """
    Build the summary cube with rankit expression sum, nnz (num cells with non zero expression) values for
    each gene for each possible group of cell attributes (cube row)

    Returns
    -------
    cube_index: pd.DataFrame
        The index of the summary cube, with one row per possible group of cell attributes

    cube_sum: np.ndarray
        The sum of expression values for each gene for each group of cell attributes

    cube_nnz: np.ndarray
        The number of cells with non zero expression for each gene for each group of cell attributes

    cube_sqsum: np.ndarray
        The squared sum of expression values for each gene for each group of cell attributes
    """

    cell_labels, cube_index = make_cube_index(tdb_group=corpus_path, cube_dims=cube_dims)
    n_groups = len(cube_index)
    n_genes = len(gene_ontology_term_ids)

    cube_sum = np.zeros((n_groups, n_genes), dtype=np.float32)
    cube_nnz = np.zeros((n_groups, n_genes), dtype=np.uint64)
    cube_sqsum = np.zeros((n_groups, n_genes), dtype=np.float32)

    reduce_X(
        tdb_group=corpus_path,
        cube_indices=cell_labels.cube_idx.values,
        cube_sum=cube_sum,
        cube_nnz=cube_nnz,
        cube_sqsum=cube_sqsum,
    )
    return TransformResult(cube_index, cube_sum, cube_nnz, cube_sqsum)


@log_func_runtime
def reduce_X(
    *, tdb_group: str, cube_indices: np.ndarray, cube_sum: np.ndarray, cube_nnz: np.ndarray, cube_sqsum: np.ndarray
) -> None:
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
                rankit_values=result["rankit"].values,
                obs_idxs=result["obs_idx"].values,
                var_idx=result["var_idx"].values,
                cube_indices=cube_indices,
                sum_into=cube_sum,
                nnz_into=cube_nnz,
                sqsum_into=cube_sqsum,
            )


# TODO: this could be further optimize by parallel chunking.  Might help large arrays if compute ends up being a bottleneck. # noqa E501
@nb.njit(fastmath=True, error_model="numpy", parallel=False, nogil=True)
def gene_expression_sum_x_cube_dimension(
    *,
    rankit_values: np.ndarray,
    obs_idxs: np.ndarray,
    var_idx: np.ndarray,
    cube_indices: np.ndarray,
    sum_into: np.ndarray,
    nnz_into: np.ndarray,
    sqsum_into: np.ndarray,
) -> None:
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
            nnz_into[grp_idx, cidx] += 1
            sqsum_into[grp_idx, cidx] += val**2


def make_cube_index(*, tdb_group: str, cube_dims: list) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create index for queryable dimensions
    """
    cell_labels = extract_obs_data(tdb_group, cube_dims)
    # number of cells with specific tuple of dims
    cube_index = pd.DataFrame(cell_labels.value_counts(), columns=["n"])

    # filter out any groups with filter_cells=True
    cube_index = cube_index[np.logical_not(list(cube_index.index.get_level_values("filter_cells")))]

    # drop filter_cells
    cube_index.index = cube_index.index.droplevel("filter_cells")

    # add cube_idx column
    cube_index["cube_idx"] = range(len(cube_index))
    cube_index["cube_idx"] = cube_index["cube_idx"].astype("int")

    # join cube_idx to cell_labels
    cell_labels = cell_labels.join(cube_index.cube_idx, on=cube_dims)

    # filter out any cells with no cube_idx, these were filtered out
    cell_labels = cell_labels[~np.isnan(cell_labels.cube_idx)]
    cell_labels["cube_idx"] = cell_labels["cube_idx"].astype("int")

    return cell_labels, cube_index
