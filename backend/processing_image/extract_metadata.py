import sys

import numpy as np
import scanpy

adata = scanpy.read_h5ad(sys.argv[1], backed="r")

print(f"organism: {adata.uns['organism']}")
print(f"tissue: {list(adata.obs.tissue.unique())}")
print(f"assay: {list(adata.obs.assay.unique())}")
print(f"disease: {list(adata.obs.disease.unique())}")
print(f"sex: {list(adata.obs.sex.unique())}")
print(f"ethnicity: {list(adata.obs.ethnicity.unique())}")
print(f"development_stage: {list(adata.obs.development_stage.unique())}")

try:
    raw_layer_name = [k for k, v in adata.uns["layer_descriptions"].items() if v == "raw"][0]
except (KeyError, IndexError):
    raise RuntimeError("Raw layer not found in layer descriptions!")

if raw_layer_name == "X":
    raw_layer = adata.X
elif raw_layer_name == "raw.X":
    raw_layer = adata.raw.X
else:
    raw_layer = adata.layers[raw_layer]

# Calling np.count_nonzero on and h5py.Dataset appears to read the entire thing
# into memory, so we need to chunk it to be safe.
stride = 25000
numerator, denominator = 0, 0
for bounds in zip(range(0, raw_layer.shape[0], stride), range(stride, raw_layer.shape[0] + stride, stride)):
    chunk = raw_layer[bounds[0] : bounds[1], :]
    numerator += np.count_nonzero(chunk)
    denominator += chunk.shape[0]

print(f"mean_genes_per_cell: {numerator/denominator}")
