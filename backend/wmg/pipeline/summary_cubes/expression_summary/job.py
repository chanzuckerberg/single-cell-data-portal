import logging

import numpy as np
import pandas as pd
import tiledb

from backend.wmg.data.schemas.cube_schema import (
    expression_summary_indexed_dims_no_gene_ontology,
    expression_summary_non_indexed_dims,
    expression_summary_schema,
)
from backend.wmg.data.schemas.cube_schema_default import (
    expression_summary_indexed_dims_no_gene_ontology as expression_summary_indexed_dims_no_gene_ontology_default,
)
from backend.wmg.data.schemas.cube_schema_default import (
    expression_summary_non_indexed_dims as expression_summary_non_indexed_dims_default,
)
from backend.wmg.data.schemas.cube_schema_default import expression_summary_schema as expression_summary_default_schema
from backend.wmg.data.snapshot import EXPRESSION_SUMMARY_CUBE_NAME, EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.data.utils import create_empty_cube, log_func_runtime
from backend.wmg.pipeline.summary_cubes.expression_summary.load import build_in_mem_cube, build_in_mem_cube_default
from backend.wmg.pipeline.summary_cubes.expression_summary.transform import transform
from backend.wmg.pipeline.summary_cubes.extract import extract_var_data

logger = logging.getLogger(__name__)


def _load(
    *,
    uri: str,
    gene_ontology_term_ids: list,
    cube_index: pd.DataFrame,
    cube_sum: np.ndarray,
    cube_nnz: np.ndarray,
    cube_sqsum: np.ndarray,
    default: bool = False,
) -> None:
    """
    Build expression summary cube in memory and write to disk
    """
    if default:
        non_indexed_dims = expression_summary_non_indexed_dims_default
        dims, vals = build_in_mem_cube_default(
            gene_ids=gene_ontology_term_ids,
            cube_index=cube_index,
            other_cube_attrs=non_indexed_dims,
            cube_sum=cube_sum,
            cube_nnz=cube_nnz,
            cube_sqsum=cube_sqsum,
        )
    else:
        non_indexed_dims = expression_summary_non_indexed_dims
        dims, vals = build_in_mem_cube(
            gene_ids=gene_ontology_term_ids,
            cube_index=cube_index,
            other_cube_attrs=non_indexed_dims,
            cube_sum=cube_sum,
            cube_nnz=cube_nnz,
        )

    logger.debug("Saving cube to tiledb")
    with tiledb.open(uri, "w") as cube:
        cube[tuple(dims)] = vals

    logger.debug("Cube created, start consolidation")
    tiledb.consolidate(uri)

    logger.debug("Cube consolidated, start vacuumming")
    tiledb.vacuum(uri)


@log_func_runtime
def create_expression_summary_cube(corpus_path: str, default=False) -> None:
    """
    Create queryable cube and write to disk
    """
    if default:
        uri = f"{corpus_path}/{EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME}"
        cube_dims = (
            expression_summary_indexed_dims_no_gene_ontology_default + expression_summary_non_indexed_dims_default
        )
    else:
        uri = f"{corpus_path}/{EXPRESSION_SUMMARY_CUBE_NAME}"
        cube_dims = expression_summary_indexed_dims_no_gene_ontology + expression_summary_non_indexed_dims
        cube_dims = [dim for dim in cube_dims if dim != "publication_citation"]

    ctx = create_ctx()

    with tiledb.scope_ctx(ctx):
        # Create cube
        if default:
            create_empty_cube(uri, expression_summary_default_schema)
        else:
            create_empty_cube(uri, expression_summary_schema)

        # extract data
        gene_ontology_term_ids = extract_var_data(corpus_path, ctx)

        # transform
        result = transform(corpus_path=corpus_path, gene_ontology_term_ids=gene_ontology_term_ids, cube_dims=cube_dims)
        _load(
            uri=uri,
            gene_ontology_term_ids=gene_ontology_term_ids,
            cube_index=result.cube_index,
            cube_sum=result.cube_sum,
            cube_nnz=result.cube_nnz,
            cube_sqsum=result.cube_sqsum,
            default=default,
        )
    gene_count = len(gene_ontology_term_ids)
    logger.info(f"create_expression_summary_cube: gene_count={gene_count}")
