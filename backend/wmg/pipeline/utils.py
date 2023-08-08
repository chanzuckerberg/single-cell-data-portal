import gc
import logging

import numpy as np
from anndata import AnnData
from scipy.sparse import csc_matrix, csr_matrix

from backend.wmg.pipeline.errors import UnsupportedMatrixTypeError

logger = logging.getLogger(__name__)

####################################### PUBLIC FUNCTIONS #########################################
def remap_anndata_normalized_X_to_raw_X_if_exists(adata: AnnData) -> None:
    """
    Remaps AnnData.X to AnnData.raw.X if and only if AnnData.raw.X exists.

    Furthermore, garbage collect the original value pointed to by AnnData.X.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    """
    raw_expression_matrix = getattr(adata.raw, "X", None)

    if raw_expression_matrix is not None:
        logger.info("Remapping AnnData.X to AnnData.raw.X")
        adata.X = adata.raw.X
        gc.collect()


def anndata_filter_cells_by_gene_counts_inplace(adata: AnnData, min_genes: int) -> None:
    """
    Filter cell outliers based on numbers of genes expressed.

    NOTE: This operation is *in place*

    This implementation is more memory efficient than `scanpy.pp.filter_cells`:
    https://github.com/scverse/scanpy/blob/89804c26bcbdde4449dfe98b3193a9d28c2a0e2f/scanpy/preprocessing/_simple.py#L41C2-L41C2

    The scanpy implementation is more memory intensive because it makes unnecessary copies
    of the underlying data structures when performing computations.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    min_genes
        Minimum number of genes expressed required for a cell to pass filtering.

    """

    # memory efficient implementation for getting count of non-zeros per row.
    # For CSR we use the `indptr` structure to compute number of non-zeros per row
    # For CSC we use the `indices` structure to compute number of non-zeros per row
    # And finally for ndarray, we simply filter the matrix and sum the resulting boolean matrix.
    if isinstance(adata.X, csr_matrix):
        logger.info(f"Filtering cells for anndata.X csr_matrix of shape: {adata.X.shape}")
        number_per_cell = np.diff(adata.X.indptr)
    elif isinstance(adata.X, csc_matrix):
        cell_indices, counts = np.unique(adata.X.indices, return_counts=True)

        # We take care of special case that an entire row has no values (implicit zero)
        # and therefore omitted in the unique count of `adata.X.indices` by allocating
        # a new array of size equal to the number of rows in the sparse matrix where each
        # index in the row corresponds to a row index in the sparse matrix and then filling
        # this new array with counts non-zero columns for a corresponding row index
        #
        # Example: 4X4 sparse matrix stored in column orientation where ROW 3 has no values.
        # That is, ROW 3 is not stored in the `adata.X.indices` array:
        # column_oriented_sparse_matrix = np.array([[1, 0, 0, 2], [0, 4, 1, 0], [0, 0, 0, 0], [0, 0, 5, 0]])
        logger.info(f"Filtering cells for anndata.X csc_matrix of shape: {adata.X.shape}")
        number_per_cell = np.zeros(adata.shape[0], dtype="int")
        number_per_cell[cell_indices] = counts
    elif isinstance(adata.X, np.ndarray):
        logger.info(f"Filtering cells for anndata.X ndarray of shape: {adata.X.shape}")
        number_per_cell = (adata.X > 0).sum(axis=1).flatten()
    else:
        raise UnsupportedMatrixTypeError(f"Unsupported AnnData.X matrix type: {type(adata.X)}")

    logger.info(f"Shape of number_per_cell array holding counts of non-zero column values: {number_per_cell.shape}")
    cell_subset = number_per_cell >= min_genes

    adata.obs["filter_cells"] = ~cell_subset

    # No subsetting is needed if there are no False values in the boolean index array.
    # The False values are the ones to filter out. When the boolean index contains only
    # True values, then the mean of the boolean index array is exactly 1.
    if cell_subset.mean() < 1:
        logger.info("Performing in place zeroing of AnnData.X rows to be filtered out")
        rows_to_zero = np.where(~cell_subset)[0]
        if isinstance(adata.X, csr_matrix):
            for row in rows_to_zero:
                adata.X.data[adata.X.indptr[row] : adata.X.indptr[row + 1]] = 0
            adata.X.eliminate_zeros()
        elif isinstance(adata.X, csc_matrix):
            adata.X.data[np.in1d(adata.X.indices, rows_to_zero)] = 0
            adata.X.eliminate_zeros()
        else:
            adata.X[~cell_subset] = 0

    else:
        logger.info("Skipping inplace zeroing of AnnData.X")
