# This module contains the schema components for the precomputed marker gene cube
# (marker_genes_cube).

import numpy as np
import tiledb

filters = [tiledb.ZstdFilter(level=+22)]

# These are the queryable cube dimensions that will be modeled as
# TileDB `Dim`s and thus can be used for _efficiently_ querying
# (slicing) the TileDB array. Order matters here!
marker_genes_indexed_dims = [
    "tissue_ontology_term_id",
    "cell_type_ontology_term_id",
]

# The cube attributes that comprise the core data stored within the cube.
# Unlike the other cube schemas, this cube schema does not have any
# distinction between logical and physical attributes.
marker_genes_attrs = [
    tiledb.Attr(name="gene_ontology_term_id", dtype="ascii", var=True, filters=filters),
    tiledb.Attr(name="marker_score", dtype=np.float32, filters=filters),
]

marker_genes_domain = tiledb.Domain(
    [
        tiledb.Dim(name=marker_genes_indexed_dim, domain=None, tile=None, dtype="ascii", filters=filters)
        for marker_genes_indexed_dim in marker_genes_indexed_dims
    ]
)

marker_genes_schema = tiledb.ArraySchema(
    domain=marker_genes_domain,
    sparse=True,
    allows_duplicates=True,
    attrs=marker_genes_attrs,
    cell_order="row-major",
    tile_order="row-major",
    capacity=10000,
)
