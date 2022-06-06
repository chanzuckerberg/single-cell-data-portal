import numpy as np
import scipy as sc


def rankit(Xraw: sc.sparse.spmatrix, offset: float = 3.0) -> sc.sparse.csr_matrix:
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
    X = Xraw.tocsr(copy=True)  # get Compressed Sparse Row format of raw expression values matrix
    indptr = X.indptr  # get row count
    for row in range(0, indptr.shape[0] - 1):
        data = X.data[indptr[row] : indptr[row + 1]]

        # Assign ranks to data, assigning the same value to ties
        ranks = sc.stats.rankdata(data, method="dense")

        max_rank = max(ranks)
        prob_level = []

        for i in ranks:
            prob_level.append(np.round((i - 0.5) / max_rank, 5))

        normal_quantiles = sc.stats.norm.ppf(prob_level, loc=offset)
        X.data[indptr[row] : indptr[row + 1]][ranks] = normal_quantiles
    return X
