import logging

import numpy as np
import pandas as pd
import tiledb

from backend.wmg.data.schemas.expression_summary_fmg_cube_schema import (
    expression_summary_fmg_indexed_dims,
    expression_summary_fmg_non_indexed_dims,
    expression_summary_fmg_schema,
)
from backend.wmg.data.snapshot import EXPRESSION_SUMMARY_FMG_CUBE_NAME
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.data.utils import create_empty_cube, log_func_runtime
from backend.wmg.pipeline.summary_cubes.expression_summary_fmg.load import build_in_mem_cube
from backend.wmg.pipeline.summary_cubes.expression_summary_fmg.transform import transform
from backend.wmg.pipeline.summary_cubes.extract import extract_var_data

logger: logging.Logger = logging.getLogger(__name__)


def _load(
    uri: str,
    gene_ontology_term_ids: list,
    cube_index: pd.DataFrame,
    cube_sum: np.ndarray,
    cube_sqsum: np.ndarray,
    cube_nnz: np.ndarray,
) -> None:
    """
    Build expression summary fmg cube in memory and write to disk
    """
    dims, vals = build_in_mem_cube(
        gene_ontology_term_ids,
        cube_index,
        expression_summary_fmg_non_indexed_dims,
        cube_sum,
        cube_sqsum,
        cube_nnz,
    )

    logger.debug("Saving cube to tiledb")
    with tiledb.open(uri, "w") as cube:
        cube[tuple(dims)] = vals

    logger.debug("Cube created, start consolidation")
    tiledb.consolidate(uri)

    logger.debug("Cube consolidated, start vacuumming")
    tiledb.vacuum(uri)


@log_func_runtime
def create_expression_summary_fmg_cube(corpus_path: str):
    """
    Create queryable cube and write to disk
    """
    uri = f"{corpus_path}/{EXPRESSION_SUMMARY_FMG_CUBE_NAME}"
    ctx = create_ctx()
    cube_dims = expression_summary_fmg_indexed_dims

    with tiledb.scope_ctx(ctx):
        # Create cube
        create_empty_cube(uri, expression_summary_fmg_schema)

        # extract data
        gene_ontology_term_ids = extract_var_data(corpus_path, ctx)

        # transform
        cube_index, cube_sum, cube_sqsum, cube_nnz = transform(corpus_path, gene_ontology_term_ids, cube_dims)
        _load(uri, gene_ontology_term_ids, cube_index, cube_sum, cube_sqsum, cube_nnz)
    gene_count = len(gene_ontology_term_ids)
    logger.info(f"create_expression_summary_fmg_cube: gene_count={gene_count}")
