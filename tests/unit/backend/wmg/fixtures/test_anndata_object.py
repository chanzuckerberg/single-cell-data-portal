import anndata
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix


def create_anndata_test_object(num_genes: int = 3, num_cells: int = 5):
    """
    Notes on csr_matrix from https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)
    csr = compressed sparse row matrix
    For example, the matrix
    [
        [5,0,0,0],
        [0,8,0,0],
        [0,0,3,0],
        [0,0,0,6]
                    ]
    is a 4 × 4 matrix with 4 nonzero elements, hence
    V         = [ 5 8 3 6 ]
    COL_INDEX = [ 0 1 2 1 ]
    ROW_INDEX = [ 0 1 2 3 4 ]
    But then there are further ways to optimize including using the indeptr
    https://stackoverflow.com/questions/52299420/scipy-csr-matrix-understand-indptr
    From scipy docs on csr_matrix
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
    Sparse matrices can be used in arithmetic operations: they support addition, subtraction,
    multiplication, division, and matrix power.
    Advantages of the CSR format
    - efficient arithmetic operations CSR + CSR, CSR * CSR, etc.
    - efficient row slicing
    - fast matrix vector products
    Disadvantages of the CSR format
    - slow column slicing operations (consider CSC)
    - changes to the sparsity structure are expensive (consider LIL or DOK)
    """
    # make matrix of gene expression (one value per gene/cell combo)
    counts = csr_matrix(np.random.poisson(1, size=(num_cells, num_genes)), dtype=np.float32)
    adata = ad.AnnData(counts)
    adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
    adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]

    # add cell level meta data
    ct = np.random.choice(["B", "T", "Monocyte"], size=(adata.n_obs,))
    adata.obs["cell_type"] = pd.Categorical(ct)  # Categoricals are preferred for efficiency
    adata.obs

    # Add cell level metadata matrices
    # matrix umap embeding
    adata.obsm["X_umap"] = np.random.normal(0, 1, size=(adata.n_obs, 2))
    # gene level metadata matrices
    adata.varm["gene_stuff"] = np.random.normal(0, 1, size=(adata.n_vars, 5))

    # unstructured metadata attached to the whole, schema version required
    adata.uns["schema_version"] = "2.0.0"

    # Add layers
    adata.layers["log_transformed"] = np.log1p(adata.X)

    return adata
