import sys
import scanpy as sc

adata = sc.read_h5ad(sys.argv[1])

# AnnData permits slashes in column names, but loom does not
column_name_map = {}
for column in adata.obs.columns:
    if "/" in column:
        column_name_map[column] = column.replace("/", "-")
if column_name_map:
    adata.obs = adata.obs.rename(columns=column_name_map)
adata.write_loom(sys.argv[1].replace(".h5ad", ".loom"), True)
