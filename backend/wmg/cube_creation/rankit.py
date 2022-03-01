import numpy as np
import scipy.stats, scipy.sparse
import numba as nb


@nb.jit
def quantiles(n: int) -> np.ndarray:
    return np.array([np.round((i - 0.5) / n, 5) for i in range(1, n + 1)])


def rankit(Xraw: scipy.sparse.spmatrix) -> scipy.sparse.csr_matrix:
    """
    https://en.wikipedia.org/wiki/Rankit

    Caveat: equal values are ranked in undefined order.
    """
    X = Xraw.tocsr(copy=True)
    indptr = X.indptr
    for row in range(0, indptr.shape[0] - 1):
        data = X.data[indptr[row] : indptr[row + 1]]
        normal_quantiles = scipy.stats.norm.ppf(quantiles(len(data)))
        rank = np.argsort(data)
        X.data[indptr[row] : indptr[row + 1]][rank] = normal_quantiles

    return X
