import numpy as np
import tiledb

from backend.wmg.data.schemas.tiledb_filters import filters_categorical, filters_numeric

expression_summary_indexed_dims = [
    "group_id",
]


expression_summary_non_indexed_dims = [
    "gene_ontology_term_id",
]

# The full set of logical cube dimensions by which the cube can be queried.
expression_summary_logical_dims = expression_summary_indexed_dims + expression_summary_non_indexed_dims


expression_summary_domain = tiledb.Domain(
    [
        tiledb.Dim(name=cube_indexed_dim, domain=None, tile=None, dtype="int", filters=filters_numeric)
        for cube_indexed_dim in expression_summary_indexed_dims
    ]
)

# The cube attributes that comprise the core data stored within the cube.
expression_summary_logical_attrs = [
    tiledb.Attr(name="sum", dtype=np.float32, filters=filters_numeric),
    tiledb.Attr(name="sqsum", dtype=np.float32, filters=filters_numeric),
]

# The TileDB `Attr`s of the cube TileDB Array. This includes the
# logical cube attributes, above, along with the non-indexed logical
# cube dimensions, which we models as TileDB `Attr`s.
expression_summary_physical_attrs = [
    tiledb.Attr(name=nonindexed_dim, dtype="ascii", var=True, filters=filters_categorical)
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
    "cell_type_ontology_term_id",
    "tissue_ontology_term_id",
    "organism_ontology_term_id",
]

cell_counts_non_indexed_dims = [
    "publication_citation",
    "disease_ontology_term_id",
    "self_reported_ethnicity_ontology_term_id",
    "sex_ontology_term_id",
]

cell_counts_logical_dims = cell_counts_indexed_dims + cell_counts_non_indexed_dims

cell_counts_domain = tiledb.Domain(
    [
        tiledb.Dim(name=cell_counts_indexed_dim, domain=None, tile=None, dtype="ascii", filters=filters_categorical)
        for cell_counts_indexed_dim in cell_counts_indexed_dims
    ]
)

cell_counts_logical_attrs = [
    # total count of cells, regardless of expression level
    tiledb.Attr(name="n_cells", dtype=np.uint32, filters=filters_numeric),
    tiledb.Attr(name="group_id", dtype=np.uint32, filters=filters_numeric),
]

cell_counts_physical_attrs = [
    tiledb.Attr(name=nonindexed_dim, dtype="ascii", var=True, filters=filters_categorical)
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
