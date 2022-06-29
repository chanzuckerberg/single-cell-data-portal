import logging

import tiledb

from backend.corpus_asset_pipelines.expression_summary_cube.extract import extract_var_data
from backend.corpus_asset_pipelines.expression_summary_cube.load import build_in_mem_cube
from backend.corpus_asset_pipelines.expression_summary_cube.transform import transform
from backend.wmg.data.schemas.cube_schema import expression_summary_schema
from backend.wmg.data.snapshot import EXPRESSION_SUMMARY_CUBE_NAME
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.data.utils import log_func_runtime, create_empty_cube
from backend.wmg.data.schemas.cube_schema import cube_non_indexed_dims, cube_indexed_dims_no_gene_ontology

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def _load(uri, gene_ontology_term_ids, cube_index, cube_sum, cube_nnz):
    dims, vals = build_in_mem_cube(gene_ontology_term_ids, cube_index, cube_non_indexed_dims, cube_sum, cube_nnz)

    logger.debug("Saving cube to tiledb")
    with tiledb.open(uri, "w") as cube:
        cube[tuple(dims)] = vals

    logger.debug("Cube created, start consolidation")
    tiledb.consolidate(uri)

    logger.debug("Cube consolidated, start vacuumming")
    tiledb.vacuum(uri)


@log_func_runtime
def create_expression_summary_cube(corpus_path):
    """
    Create queryable cube and write to disk
    """
    uri = f"{corpus_path}/{EXPRESSION_SUMMARY_CUBE_NAME}"
    ctx = create_ctx()
    cube_dims = cube_indexed_dims_no_gene_ontology + cube_non_indexed_dims

    with tiledb.scope_ctx(ctx):
        # Create cube
        create_empty_cube(uri, expression_summary_schema)

        # extract data
        gene_ontology_term_ids = extract_var_data(corpus_path, ctx)

        # transform
        cube_index, cube_sum, cube_nnz = transform(corpus_path, gene_ontology_term_ids, cube_dims)

        _load(uri, gene_ontology_term_ids, cube_index, cube_sum, cube_nnz)
