import numpy as np
import tiledb

filters = [tiledb.ZstdFilter(level=+22)]

# Cell Counts Array

cell_counts_indexed_dims = [
    "cell_type_ontology_term_id",
    "tissue_ontology_term_id",
    "organism_ontology_term_id",
]
cell_counts_non_indexed_dims = [
    "dataset_id",
    "disease_ontology_term_id",
    "self_reported_ethnicity_ontology_term_id",
    "sex_ontology_term_id",
]


cell_counts_logical_dims = cell_counts_indexed_dims + cell_counts_non_indexed_dims


cell_counts_domain = tiledb.Domain(
    [
        tiledb.Dim(name=cell_counts_indexed_dim, domain=None, tile=None, dtype="ascii", filters=filters)
        for cell_counts_indexed_dim in cell_counts_indexed_dims
    ]
)

cell_counts_logical_attrs = [
    tiledb.Attr(name="n_cells", dtype=np.uint32, filters=filters),
]

cell_counts_physical_attrs = [
    tiledb.Attr(name=nonindexed_dim, dtype="ascii", var=True, filters=filters)
    for nonindexed_dim in cell_counts_non_indexed_dims
] + cell_counts_logical_attrs


cell_counts_schema = tiledb.ArraySchema(
    domain=cell_counts_domain,
    sparse=True,
    allows_duplicates=True,
    attrs=cell_counts_physical_attrs,
    cell_order="row-major",
    tile_order="row-major",
    capacity=10000,
)
