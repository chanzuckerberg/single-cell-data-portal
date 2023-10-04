import logging

import numba as nb
import numpy as np
from scipy import sparse, stats

logger = logging.getLogger("wmg")


@nb.jit
def quantiles(max_rank: int, ranks: np.ndarray) -> np.ndarray:
    """
    :returns an array of n floats equally spaced from 0 to 1
    """
    return np.array([np.round((i - 0.5) / max_rank, 5) for i in ranks])


def rankit(X: sparse.csr_matrix, offset: float = 3.0):
    """
    Row-wise normalizes values of a matrix using the rankit method. The target distribution is a normal distribution
    with variance of 1 and mean as set in `offset`
    https://en.wikipedia.org/wiki/Rankit
    In statistics, rankits of a set of data are the expected values of the order statistics of
    a sample from the standard normal distribution the same size as the data
    Caveat: equal values are ranked in undefined order.
    param Xraw: query matrix to be normalized
    param offset: mean for the resulting row-wise values that will follow a normal distribution with variance 1. This
    helps to shift values to a positive scale.
    :returns row-wise normalized matrix using rankit
    """
    indptr = X.indptr  # get row count
    warning_raised = False
    for row in range(0, indptr.shape[0] - 1):
        data = X.data[indptr[row] : indptr[row + 1]]
        if len(data) > 0:
            # Assign ranks to data, assigning the same value to ties
            ranks = stats.rankdata(data, method="dense")

            max_rank = max(ranks)
            prob_level = quantiles(max_rank, ranks)

            normal_quantiles = stats.norm.ppf(prob_level, loc=offset)
            X.data[indptr[row] : indptr[row + 1]] = normal_quantiles
        elif not warning_raised:
            logging.warn("This dataset has at least one row of all zero expressions")
            warning_raised = True
