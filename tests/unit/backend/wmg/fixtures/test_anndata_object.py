import os
import pathlib

import anndata
import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from backend.wmg.data.constants import INCLUDED_ASSAYS


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
    # create sparse matrix of gene expression (one value per gene/cell combo)
    # values are randomly populated using a poisson distribution
    counts = csr_matrix(np.random.poisson(1, size=(num_cells, num_genes)), dtype=np.float32)
    adata = ad.AnnData(counts)
    adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
    adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]

    # add cell level metadata
    cell_type = np.random.choice(["B", "T", "Monocyte"], size=(adata.n_obs,))
    adata.obs["cell_type"] = pd.Categorical(cell_type)  # Categoricals are preferred for efficiency
    assay_ontologies = np.random.choice(list(INCLUDED_ASSAYS.keys()), size=(adata.n_obs,))
    adata.obs["assay_ontology_term_id"] = pd.Categorical(assay_ontologies)
    adata.obs["tissue_ontology_term_id"] = "UBERON:0000101"
    adata.obs["tissue"] = "lobe of lung"
    adata.obs["tissue_ontology_term_id"] = pd.Categorical(adata.obs["tissue_ontology_term_id"])
    adata.obs["tissue"] = pd.Categorical(adata.obs["tissue"])

    # Add cell level metadata matrices
    # matrix umap embeding
    adata.obsm["X_umap"] = np.random.normal(0, 1, size=(adata.n_obs, 2))
    # gene level metadata matrices
    adata.varm["gene_stuff"] = np.random.normal(0, 1, size=(adata.n_vars, 5))

    # unstructured metadata attached to the whole, schema version required
    adata.uns["schema_version"] = "3.0.0"

    # Add layers
    adata.layers["log_transformed"] = np.log1p(adata.X)

    return adata


def store_anndata_test_file(path, dataset_name, annadata_object):
    dataset_path = f"{path}/{dataset_name}"
    os.mkdir(dataset_path)
    test_anndata_file_name = pathlib.Path(dataset_path, "local.h5ad")
    test_anndata_file_name.touch()
    annadata_object.write(test_anndata_file_name, compression="gzip")
    return test_anndata_file_name


def create_anndata_test_fixture(path: str, dataset_name: str, num_genes: int = 1000, num_cells: int = 5000) -> str:
    test_anndata_object = create_anndata_test_object(num_genes=num_genes, num_cells=num_cells)
    test_anndata_file_name = store_anndata_test_file(path, dataset_name, test_anndata_object)
    return test_anndata_file_name
