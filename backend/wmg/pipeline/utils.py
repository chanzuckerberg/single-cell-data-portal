import logging
import operator
from typing import Optional, Tuple

import numpy as np
from anndata import AnnData
from scanpy.preprocessing._distributed import materialize_as_ndarray
from scipy.sparse import issparse

logger = logging.getLogger(__name__)

####################################### PUBLIC FUNCTIONS #########################################
def anndata_filter_cells(
    data: AnnData,
    min_counts: Optional[int] = None,
    min_genes: Optional[int] = None,
    max_counts: Optional[int] = None,
    max_genes: Optional[int] = None,
    inplace: bool = True,
    copy: bool = False,
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """
    Filter cell outliers based on counts and numbers of genes expressed.

    For instance, only keep cells with at least `min_counts` counts or
    `min_genes` genes expressed. This is to filter measurement outliers,
    i.e. “unreliable” observations.

    Only provide one of the optional parameters `min_counts`, `min_genes`,
    `max_counts`, `max_genes` per call.

    ###################### IMPORTANT NOTE ################################

    This is a patched version of `scanpy.pp.filter_cells`:
    https://github.com/scverse/scanpy/blob/89804c26bcbdde4449dfe98b3193a9d28c2a0e2f/scanpy/preprocessing/_simple.py#L41C2-L41C2

    The scanpy implementation is more memory intensive due to copies
    made of the underyling structures supporting a compressed sparse matrix
    when performing a filtering operation on the matrix by comparing a matrix
    with a scalar value. Particularly, the following line of code allocates additional
    memory when performing the comparison AnnData.X > 0:

    https://github.com/scverse/scanpy/blob/89804c26bcbdde4449dfe98b3193a9d28c2a0e2f/scanpy/preprocessing/_simple.py#L146

    This comparison is being replaced by the implementation in the private function:
    _filter_compressed_matrix_gt_zero()

    #######################################################################

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    min_counts
        Minimum number of counts required for a cell to pass filtering.
    min_genes
        Minimum number of genes expressed required for a cell to pass filtering.
    max_counts
        Maximum number of counts required for a cell to pass filtering.
    max_genes
        Maximum number of genes expressed required for a cell to pass filtering.
    inplace
        Perform computation inplace or return result.

    Returns
    -------
    Depending on `inplace`, returns the following arrays or directly subsets
    and annotates the data matrix:

    cells_subset
        Boolean index mask that does filtering. `True` means that the
        cell is kept. `False` means the cell is removed.
    number_per_cell
        Depending on what was thresholded (`counts` or `genes`),
        the array stores `n_counts` or `n_cells` per gene.
    """
    if copy:
        logger.warning("`copy` is deprecated, use `inplace` instead.")
    n_given_options = sum(option is not None for option in [min_genes, min_counts, max_genes, max_counts])
    if n_given_options != 1:
        raise ValueError(
            "Only provide one of the optional parameters `min_counts`, "
            "`min_genes`, `max_counts`, `max_genes` per call."
        )
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        cell_subset, number = materialize_as_ndarray(
            anndata_filter_cells(adata.X, min_counts, min_genes, max_counts, max_genes)
        )
        if not inplace:
            return cell_subset, number
        if min_genes is None and max_genes is None:
            adata.obs["n_counts"] = number
        else:
            adata.obs["n_genes"] = number
        adata._inplace_subset_obs(cell_subset)
        return adata if copy else None
    X = data  # proceed with processing the data matrix
    min_number = min_counts if min_genes is None else min_genes
    max_number = max_counts if max_genes is None else max_genes

    ################ IMPORTANT NOTE ##########################
    # This the implementation change from `scanpy.pp.filter_cells()`.
    # That is the following line is replaced by a custom implemention
    # of filtering the compressed matrix
    # https://github.com/scverse/scanpy/blob/89804c26bcbdde4449dfe98b3193a9d28c2a0e2f/scanpy/preprocessing/_simple.py#L146
    filtered_matrix = X
    if not (min_genes is None and max_genes is None):
        filtered_matrix = _filter_compressed_matrix_gt_zero(X)

    number_per_cell = np.sum(filtered_matrix, axis=1)
    #############################################################
    if issparse(X):
        number_per_cell = number_per_cell.A1
    if min_number is not None:
        cell_subset = number_per_cell >= min_number
    if max_number is not None:
        cell_subset = number_per_cell <= max_number

    s = np.sum(~cell_subset)
    if s > 0:
        msg = f"filtered out {s} cells that have "
        if min_genes is not None or min_counts is not None:
            msg += "less than "
            msg += f"{min_genes} genes expressed" if min_counts is None else f"{min_counts} counts"
        if max_genes is not None or max_counts is not None:
            msg += "more than "
            msg += f"{max_genes} genes expressed" if max_counts is None else f"{max_counts} counts"
        logger.info(msg)
    return cell_subset, number_per_cell


################################## PRIVATE FUNCTIONS ############################################


def _filter_compressed_matrix_gt_zero(csr_matrix):
    """
    Filters a compressed matrix such that the filtered matrix has elements > 0.

    The implementation is almost the same as `scipy.sparse.compressed._scalar_binopt()`,
    except with modifications to not make copies of structure arrays underlying the matrix:

    https://github.com/scipy/scipy/blob/c025587b9fbbabf7e6bbf10e2e7c76ae61d6395c/scipy/sparse/_compressed.py#L208

    Parameters
    ----------
    csr_matrix: This is an instance of `scipy.sparse.csr_matrix` class
    """
    csr_matrix.sum_duplicates()

    # Setting `copy=False` is the key implementation change from
    # `scipy.sparse.compressed._scalar_binopt()`. The scipy implementation
    # sets `copy=True` which creates copies. We avoid creating copies here.
    res = csr_matrix._with_data(operator.gt(csr_matrix.data, 0), copy=False)
    res.eliminate_zeros()
    return res
