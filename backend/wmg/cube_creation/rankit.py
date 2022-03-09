import numpy as np
import scipy.stats, scipy.sparse
import numba as nb


@nb.jit
def quantiles(n: int) -> np.ndarray:
    """
    :returns an array of n floats equally spaced from 0 to 1
    """
    return np.array(
        [np.round((i - 0.5) / n, 5) for i in range(1, n + 1)]
    )


def rankit(Xraw: scipy.sparse.spmatrix, offset: float = 3.0) -> scipy.sparse.csr_matrix:
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
    X = Xraw.tocsr(copy=True) # get Compressed Sparse Row format of raw expression values matrix
    indptr = X.indptr # get row count
    for row in range(0, indptr.shape[0] - 1):
        data = X.data[indptr[row] : indptr[row + 1]]
        # A normal continuous random variable.
        normal_quantiles = scipy.stats.norm.ppf(quantiles(len(data)), loc=offset)
        rank = np.argsort(data)
        X.data[indptr[row] : indptr[row + 1]][rank] = normal_quantiles

    return X
