import os
from typing import Tuple

import numpy as np
from numba import njit, prange
from scipy import stats

from backend.common.utils.rollup import are_cell_types_colinear
from backend.wmg.data.utils import setup_retry_session


@njit(parallel=True)
def nanpercentile_2d(arr: np.ndarray, percentile: float, axis: int) -> np.ndarray:
    """
    Calculate the specified percentile of a 2D array along an axis, ignoring NaN values.

    Arguments
    ---------
    arr - 2D array to calculate percentile of
    percentile - percentile to calculate, as a number between 0 and 100
    axis - axis along which to calculate percentile

    Returns
    -------
    The specified percentile of the 2D array along the specified axis.
    """
    if axis == 0:
        result = np.empty(arr.shape[1])
        for i in prange(arr.shape[1]):
            arr_column = arr[:, i]
            result[i] = nanpercentile(arr_column, percentile)
        return result
    else:
        result = np.empty(arr.shape[0])
        for i in prange(arr.shape[0]):
            arr_row = arr[i, :]
            result[i] = nanpercentile(arr_row, percentile)
        return result


@njit
def nanpercentile(arr: np.ndarray, percentile: float):
    """
    Calculate the specified percentile of an array, ignoring NaN values.

    Arguments
    ---------
    arr - array to calculate percentile of
    percentile - percentile to calculate, as a number between 0 and 100

    Returns
    -------
    The specified percentile of the array.
    """

    arr_without_nan = arr[np.logical_not(np.isnan(arr))]
    length = len(arr_without_nan)

    if length == 0:
        return np.nan

    return np.percentile(arr_without_nan, percentile)


def run_ttest(
    *, sum1: np.ndarray, sumsq1: np.ndarray, n1: np.ndarray, sum2: np.ndarray, sumsq2: np.ndarray, n2: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Run a t-test on two sets of data, element-wise.
    Arrays "1" and "2" have to be broadcastable into each other.

    Arguments
    ---------
    sum1 - Sum of the first set of data
    sumsq1 - Sum of squares of the first set of data
    n1 - Number of elements in the first set of data
    sum2 - Sum of the second set of data
    sumsq2 - Sum of squares of the second set of data
    n2 - Number of elements in the second set of data

    Returns
    -------
    pvals - The p-values of the t-test
    effects - The effect sizes of the t-test
    """

    with np.errstate(divide="ignore", invalid="ignore"):
        mean1 = sum1 / n1
        meansq1 = sumsq1 / n1

        mean2 = sum2 / n2
        meansq2 = sumsq2 / n2

        var1 = meansq1 - mean1**2
        var1[var1 < 0] = 0
        var2 = meansq2 - mean2**2
        var2[var2 < 0] = 0

        var1_n = var1 / n1
        var2_n = var2 / n2
        sum_var_n = var1_n + var2_n
        dof = sum_var_n**2 / (var1_n**2 / (n1 - 1) + var2_n**2 / (n2 - 1))
        tscores = (mean1 - mean2) / np.sqrt(sum_var_n)
        effects = (mean1 - mean2) / np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 1))

    pvals = stats.t.sf(tscores, dof)
    return pvals, effects


def post_process_stats(
    *,
    cell_type_target: str,
    cell_types_context: np.ndarray,
    genes: np.ndarray,
    pvals: np.ndarray,
    effects: np.ndarray,
    percentile: float = 0.05,
) -> dict[str, dict[str, float]]:
    """
    Post-process the statistical results to handle colinearity of cell types in the ontology and calculate percentiles.

    Arguments
    ---------
    cell_type_target - The target cell type
    cell_types_context - The context cell types
    genes - The genes involved in the analysis
    pvals - The p-values from the statistical test
    effects - The effect sizes from the statistical test
    percentile - The percentile to use for thresholding (default is 0.05)

    Returns
    -------
    A dictionary mapping marker genes to their statistics.
    """

    is_colinear = np.array([are_cell_types_colinear(cell_type, cell_type_target) for cell_type in cell_types_context])
    effects[is_colinear] = np.nan
    pvals[is_colinear] = np.nan

    pvals[:, np.all(np.isnan(pvals), axis=0)] = 1
    effects[:, np.all(np.isnan(effects), axis=0)] = 0

    # aggregate
    effects = nanpercentile_2d(effects, percentile * 100, 0)

    effects[effects == 0] = np.nan

    # pvals = np.array([stats.combine_pvalues(x[np.invert(np.isnan(x))] + 1e-300)[-1] for x in pvals.T])
    pvals = np.sort(pvals, axis=0)[int(np.round(0.05 * pvals.shape[0]))]

    markers = np.array(genes)[np.argsort(-effects)]
    p = pvals[np.argsort(-effects)]
    effects = effects[np.argsort(-effects)]

    statistics = []
    final_markers = []
    for i in range(len(p)):
        pi = p[i]
        ei = effects[i]
        if ei is not np.nan and pi is not np.nan:
            statistics.append({"p_value": pi, "effect_size": ei})
            final_markers.append(markers[i])
    return dict(zip(list(final_markers), statistics))


def query_gene_info_for_gene_description(gene_id: str) -> str:
    """
    Query the gene description for a given gene ID.

    Arguments
    ---------
    gene_id - The ID of the gene to query

    Returns
    -------
    The name of the gene if the query is successful, otherwise returns the gene ID.
    """

    API_URL = os.getenv("API_URL")
    if API_URL is None:
        API_URL = "https://api.cellxgene.dev.single-cell.czi.technology"

    session = setup_retry_session()
    response = session.get(f"{API_URL}/gene_info/v1/gene_info?gene={gene_id}")
    if response.status_code == 200:
        data = response.json()
        return data["name"]
    else:
        return gene_id
