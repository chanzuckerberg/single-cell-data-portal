import numba as nb
import numpy as np


"""
TODO: this could be further optimize by parallel chunking.  Might help large
arrays if compute ends up being a bottlneck.
"""


@nb.njit(fastmath=True, error_model="numpy", parallel=True, nogil=True)
def mean_var_coo(n_rows, n_cols, data, col):
    mean = np.zeros((n_cols,), dtype=np.float64)
    var = np.zeros((n_cols,), dtype=np.float64)
    nonfinite_count = np.zeros((n_cols,), dtype=np.int32)
    nonzero_count = np.zeros((n_cols,), dtype=np.int32)

    # sum by axis=0
    for k in nb.prange(0, len(data)):
        val = data[k]
        cidx = col[k]
        if np.isfinite(val):
            mean[cidx] += val
        else:
            nonfinite_count[cidx] += 1

    mean /= n_rows - nonfinite_count

    for k in nb.prange(0, len(data)):
        val = data[k]
        cidx = col[k]
        if np.isfinite(val):  # ignore NaN and +/-Inf
            dfm = val - mean[cidx]
            var[cidx] += dfm * dfm
            nonzero_count[cidx] += 1

    var += (n_rows - nonzero_count) * mean * mean
    var /= n_rows - nonfinite_count - 1

    return mean, var


@nb.njit(fastmath=True, error_model="numpy", parallel=False, nogil=True)
def coo_cube_pass1_into(data, row, col, row_groups, sum_into, nnz_into, min_into, max_into):
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


@nb.njit(fastmath=True, error_model="numpy", parallel=False, nogil=True)
def coo_cube_pass2_into(data, row, col, row_groups, min_vals, max_vals, n_bins, bins_into):
    bin_widths = (max_vals - min_vals) / n_bins
    for k in range(0, len(data)):
        val = data[k]
        if np.isfinite(val):
            cidx = col[k]
            grp_idx = row_groups[row[k]]
            bin_width = bin_widths[grp_idx][cidx]
            bin = int((val - min_vals[grp_idx][cidx]) / bin_width)
            if bin >= n_bins:
                bin = n_bins - 1
            bins_into[grp_idx][cidx][bin] += 1
