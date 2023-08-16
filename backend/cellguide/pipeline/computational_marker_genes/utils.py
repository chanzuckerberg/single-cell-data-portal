import numpy as np
from numba import njit, prange
from scipy import stats

from backend.common.utils.rollup import are_cell_types_colinear


@njit(parallel=True)
def nanpercentile_2d(arr, percentile, axis):
    """
    Calculate the specified percentile of a 2D array along an axis, ignoring NaN values.

    Parameters:
        arr: 2D array to calculate percentile of
        percentile: percentile to calculate, as a number between 0 and 100
        axis: axis along which to calculate percentile

    Returns:
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
def nanpercentile(arr, percentile):
    arr_without_nan = arr[np.logical_not(np.isnan(arr))]
    length = len(arr_without_nan)

    if length == 0:
        return np.nan

    return np.percentile(arr_without_nan, percentile)


def run_ttest(sum1, sumsq1, n1, sum2, sumsq2, n2):
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


def post_process_stats(cell_type_target, cell_types_context, genes, pvals, effects, percentile=0.05):
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
