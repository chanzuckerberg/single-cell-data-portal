import anndata
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix


def create_anndata_test_object(num_genes: int = 3, num_cells: int = 5):
    # make matrix of gene expression (one value per gene/cell combo)
    counts = csr_matrix(np.random.poisson(1, size=(num_genes, num_cells)), dtype=np.float32)
    adata = ad.AnnData(counts)
    adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
    adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]

    # add cell level meta data
    ct = np.random.choice(["B", "T", "Monocyte"], size=(adata.n_obs,))
    adata.obs["cell_type"] = pd.Categorical(ct)  # Categoricals are preferred for efficiency
    adata.obs

    ## add cell level metadata matrices
    # matrix umap embeding
    adata.obsm["X_umap"] = np.random.normal(0, 1, size=(adata.n_obs, 2))
    # gene level metadata matrices
    adata.varm["gene_stuff"] = np.random.normal(0, 1, size=(adata.n_vars, 5))

    # unstructured metadata attached to the whole
    adata.uns["random"] = [1, 2, 3]

    ## add layers
    adata.layers["log_transformed"] = np.log1p(adata.X)

    return adata
