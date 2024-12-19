import logging

import dask.array as da
from cellxgene_schema.validate import Validator

logger: logging.Logger = logging.getLogger("matrix_utils")


def is_matrix_sparse(matrix: da.Array, sparse_threshold: float) -> bool:
    """
    Returns whether `matrix` is sparse or not (i.e. dense). This is determined by figuring out whether the matrix has
    a sparsity percentage below the sparse_threshold, returning the number of non-zeros encountered and number of
    elements evaluated.

    TODO move to dask_utils.py in a separate PR
    """
    if sparse_threshold == 100.0:
        return True
    if sparse_threshold == 0.0:
        return False

    total_number_of_matrix_elements = matrix.shape[0] * matrix.shape[1]
    number_of_non_zero_elements = Validator.count_matrix_nonzero(matrix)
    is_sparse = (100.0 * (number_of_non_zero_elements / total_number_of_matrix_elements)) < sparse_threshold
    return is_sparse


def enforce_canonical_format(adata):
    """
    Enforce canonical format for an AnnData, if not already in canonical format.  This function will modify the
    matrix in place.
    """
    X = adata.X
    if hasattr(X, "has_canonical_format") and not X.has_canonical_format:
        logger.warning("noncanonical data found in X; converting to canonical format using sum_duplicates.")
        X.sum_duplicates()
