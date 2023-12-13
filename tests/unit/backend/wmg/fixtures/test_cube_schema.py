import numpy as np
import tiledb

# the temporary wmg snapshot used for testing will contain data for every combination
# of the below logical dimensions. These dimensions are a subset of the dimensions
# present in the real snapshot cubes to avoid repeatedly generating unnecessarily large
# test fixtures.

expression_summary_indexed_dims = [
    "gene_ontology_term_id",
    "tissue_ontology_term_id",
    "organism_ontology_term_id",
]

expression_summary_indexed_dims_no_gene_ontology = [
    "tissue_ontology_term_id",
    "organism_ontology_term_id",
]

expression_summary_non_indexed_dims = [
    "cell_type_ontology_term_id",
    "development_stage_ontology_term_id",
    "self_reported_ethnicity_ontology_term_id",
]

expression_summary_logical_dims = expression_summary_indexed_dims + expression_summary_non_indexed_dims

filters = [tiledb.ZstdFilter(level=+22)]

expression_summary_domain = tiledb.Domain(
    [
        tiledb.Dim(name=cube_indexed_dim, domain=None, tile=None, dtype="ascii", filters=filters)
        for cube_indexed_dim in expression_summary_indexed_dims
    ]
)

expression_summary_logical_attrs = [
    tiledb.Attr(name="nnz", dtype=np.uint64, filters=filters),
    tiledb.Attr(name="sum", dtype=np.float32, filters=filters),
    tiledb.Attr(name="sqsum", dtype=np.float32, filters=filters),
]

expression_summary_physical_attrs = [
    tiledb.Attr(name=nonindexed_dim, dtype="ascii", var=True, filters=filters)
    for nonindexed_dim in expression_summary_non_indexed_dims
] + expression_summary_logical_attrs

expression_summary_schema = tiledb.ArraySchema(
    domain=expression_summary_domain,
    sparse=True,
    allows_duplicates=True,
    attrs=expression_summary_physical_attrs,
    cell_order="row-major",
    tile_order="row-major",
    capacity=10000,
)

# Cell Counts Array

cell_counts_indexed_dims = [
    "tissue_ontology_term_id",
    "organism_ontology_term_id",
]
cell_counts_non_indexed_dims = expression_summary_non_indexed_dims

cell_counts_logical_dims = cell_counts_indexed_dims + cell_counts_non_indexed_dims


cell_counts_domain = tiledb.Domain(
    [
        tiledb.Dim(name=cell_counts_indexed_dim, domain=None, tile=None, dtype="ascii", filters=filters)
        for cell_counts_indexed_dim in cell_counts_indexed_dims
    ]
)

cell_counts_logical_attrs = [
    # total count of cells, regardless of expression level
    tiledb.Attr(name="n_cells", dtype=np.uint32, filters=filters),
]

cell_counts_physical_attrs = [
    tiledb.Attr(name=nonindexed_dim, dtype="ascii", var=True, filters=filters)
    for nonindexed_dim in expression_summary_non_indexed_dims
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
