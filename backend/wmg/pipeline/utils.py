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
        number_per_cell = np.diff(adata.X.indptr)
    elif isinstance(adata.X, csc_matrix):
        _, number_per_cell = np.unique(adata.X.indices, return_counts=True)
    elif isinstance(adata.X, np.ndarray):
        number_per_cell = (adata.X > 0).sum(axis=1).A1
    else:
        raise UnsupportedMatrixTypeError(f"Unsupported AnnData.X matrix type: {type(adata.X)}")

    cell_subset = number_per_cell >= min_genes

    # No subsetting is needed if there are no False values in the boolean index array.
    # The False values are the ones to filter out. When the boolean index contains only
    # True values, then the mean of the boolean index array is exactly 1.
    if cell_subset.mean() < 1:
        adata._inplace_subset_obs(cell_subset)
