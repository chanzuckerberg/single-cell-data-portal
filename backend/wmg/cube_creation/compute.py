import numba as nb
import numpy as np


"""
TODO: this could be further optimize by parallel chunking.  Might help large
arrays if compute ends up being a bottlneck.
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
