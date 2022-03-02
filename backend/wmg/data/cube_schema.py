import numpy as np
import tiledb

cube_indexed_dims = [
    "gene_ontology_term_id",
    "tissue_ontology_term_id",
    "organism_ontology_term_id",
]

cube_non_indexed_dims = [
    "cell_type_ontology_term_id",
    "dataset_id",
    "assay_ontology_term_id",
    "development_stage_ontology_term_id",
    "disease_ontology_term_id",
    "ethnicity_ontology_term_id",
    "sex_ontology_term_id",
]
cube_logical_dims = cube_indexed_dims + cube_non_indexed_dims

filters = [tiledb.ZstdFilter(level=+22)]

domain = tiledb.Domain(
    [
        tiledb.Dim(name=cube_indexed_dim, domain=None, tile=None, dtype="ascii", filters=filters)
        for cube_indexed_dim in cube_indexed_dims
    ]
)

# summary expression data
cube_logical_attrs = [
    tiledb.Attr(name="n_cells", dtype=np.uint32, filters=filters),
    tiledb.Attr(name="nnz", dtype=np.uint64, filters=filters),  # TODO: Why uint64?
    tiledb.Attr(name="sum", dtype=np.float32, filters=filters),
]

# metadata indexes to search along
cube_physical_attrs = [
    tiledb.Attr(name=nonindexed_dim, dtype="ascii", var=True, filters=filters)
    for nonindexed_dim in cube_non_indexed_dims
] + cube_logical_attrs

schema = tiledb.ArraySchema(
    domain=domain,
    sparse=True,
    allows_duplicates=True,
    attrs=cube_physical_attrs,
    cell_order="row-major",
    tile_order="row-major",
    capacity=10000,
)
