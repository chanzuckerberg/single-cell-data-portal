# This module contains the schema components for the expression summary cube.
# used for marker gene calculation (expression_summary_fmg).

import numpy as np
import tiledb

filters = [tiledb.ZstdFilter(level=+22)]

# These are the queryable cube dimensions that will be modeled as
# TileDB `Dim`s and thus can be used for _efficiently_ querying
# (slicing) the TileDB array. Order matters here!
expression_summary_fmg_indexed_dims = [
    "tissue_ontology_term_id",
    "organism_ontology_term_id",
    "cell_type_ontology_term_id",
]

# These are the queryable cube dimensions that will be modeled as
# TileDB `Attrs` (i.e. (non-indexed") and thus will require
# client-side filtering, which may result in less efficient querying.
expression_summary_fmg_non_indexed_dims = [
    "gene_ontology_term_id",
]

# The full set of logical cube dimensions by which the cube can be queried.
expression_summary_fmg_logical_dims = expression_summary_fmg_indexed_dims + expression_summary_fmg_non_indexed_dims

expression_summary_fmg_domain = tiledb.Domain(
    [
        tiledb.Dim(name=cube_indexed_dim, domain=None, tile=None, dtype="ascii", filters=filters)
        for cube_indexed_dim in expression_summary_fmg_indexed_dims
    ]
)

# The cube attributes that comprise the core data stored within the cube.
expression_summary_fmg_logical_attrs = [
    tiledb.Attr(name="sum", dtype=np.float32, filters=filters),
    tiledb.Attr(name="sqsum", dtype=np.float32, filters=filters),
    tiledb.Attr(name="nnz", dtype=np.uint64, filters=filters),
]

# The TileDB `Attr`s of the cube TileDB Array. This includes the
# logical cube attributes, above, along with the non-indexed logical
# cube dimensions, which we models as TileDB `Attr`s.
expression_summary_fmg_physical_attrs = [
    tiledb.Attr(name=nonindexed_dim, dtype="ascii", var=True, filters=filters)
    for nonindexed_dim in expression_summary_fmg_non_indexed_dims
] + expression_summary_fmg_logical_attrs

expression_summary_fmg_schema = tiledb.ArraySchema(
    domain=expression_summary_fmg_domain,
    sparse=True,
    allows_duplicates=True,
    attrs=expression_summary_fmg_physical_attrs,
    cell_order="row-major",
    tile_order="row-major",
    capacity=10000,
)
