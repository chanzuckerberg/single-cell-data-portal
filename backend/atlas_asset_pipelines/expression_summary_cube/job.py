import logging
import time

import tiledb

from backend.atlas_asset_pipelines.expression_summary_cube.extract import extract_var_data
from backend.atlas_asset_pipelines.expression_summary_cube.transform import transform
from backend.wmg.data.schemas.cube_schema import expression_summary_schema, cube_non_indexed_dims
from backend.wmg.data.snapshot import EXPRESSION_SUMMARY_CUBE_NAME
from backend.wmg.data.tiledb import create_ctx
from backend.atlas_asset_pipelines.expression_summary_cube.load import build_in_mem_cube
from backend.wmg.data.utils import create_empty_cube

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def create_expression_summary_cube(corpus_path):
    """
    Create queryable cube and write to disk
    """
    uri = f"{corpus_path}/{EXPRESSION_SUMMARY_CUBE_NAME}"
    cube_creation_start_time = time.time()
    # todo safe te reuse ctx here?
    ctx = create_ctx()

    with tiledb.scope_ctx(ctx):
        # Create cube
        logger.info(f"Start creating expression summary cube at : {uri}")
        create_empty_cube(uri, expression_summary_schema)

        # extract data
        gene_ontology_term_ids = extract_var_data(corpus_path, ctx)

        # transform
        cube_index, cube_sum, cube_nnz = transform(corpus_path, gene_ontology_term_ids)

        # load data
        dims, vals = build_in_mem_cube(gene_ontology_term_ids, cube_index, cube_non_indexed_dims, cube_sum, cube_nnz)

        logger.debug("Saving cube to tiledb")
        with tiledb.open(uri, "w") as cube:
            cube[tuple(dims)] = vals

        logger.debug("Cube created, start consolidation")
        tiledb.consolidate(uri)

        logger.debug("Cube consolidated, start vacuumming")
        tiledb.vacuum(uri)

    logger.debug("Cube creation complete")
    create_cube_sec = time.time() - cube_creation_start_time
    logger.info(f"Expression summary cube: time to create {create_cube_sec}, uri={uri}")
