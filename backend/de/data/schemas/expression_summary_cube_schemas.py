import numpy as np
import tiledb

def generate_schema(expression_summary_indexed_dims, expression_summary_non_indexed_dims):
    filters = [tiledb.ZstdFilter(level=+22)]

    expression_summary_domain = tiledb.Domain(
        [
            tiledb.Dim(name=cube_indexed_dim, domain=None, tile=None, dtype="ascii", filters=filters)
            for cube_indexed_dim in expression_summary_indexed_dims
        ]
    )

    expression_summary_logical_attrs = [
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
        cell_order="hilbert",
        capacity=100000,
    )
    return expression_summary_schema

base_expression_summary_indexed_dims = [
    "cell_type_ontology_term_id",
    "tissue_ontology_term_id",
    "organism_ontology_term_id",
]

expression_summary_secondary_dims = [
    "dataset_id",
    "disease_ontology_term_id",
    "self_reported_ethnicity_ontology_term_id",
    "sex_ontology_term_id",
]

expression_summary_attrs = [
    "gene_ontology_term_id"
]

schemas = {"default": generate_schema(base_expression_summary_indexed_dims, expression_summary_attrs)}
for secondary_dim in expression_summary_secondary_dims:
    expression_summary_indexed_dims = base_expression_summary_indexed_dims+[secondary_dim]
    expression_summary_non_indexed_dims = [dim for dim in expression_summary_secondary_dims if dim != secondary_dim] + expression_summary_attrs
    schemas[secondary_dim] = generate_schema(expression_summary_indexed_dims, expression_summary_non_indexed_dims)    

