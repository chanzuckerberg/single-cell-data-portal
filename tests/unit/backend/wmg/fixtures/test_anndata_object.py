import anndata
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix


def create_anndata_test_object(num_genes: int = 3, num_cells: int = 5):
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
