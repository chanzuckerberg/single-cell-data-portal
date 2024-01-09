import logging
import os

import pandas as pd
import tiledb

from backend.wmg.data.schemas.cube_schema_cell_type import (
    expression_summary_schema,
)
from backend.wmg.data.snapshot import (
    EXPRESSION_SUMMARY_CELL_TYPE_CUBE_NAME,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME,
)
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_CELL_TYPE_CUBE_CREATED_FLAG,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG,
)
from backend.wmg.pipeline.errors import PipelineStepMissing
from backend.wmg.pipeline.utils import (
    create_empty_cube_if_needed,
    load_pipeline_state,
    log_func_runtime,
    write_pipeline_state,
)

logger = logging.getLogger(__name__)


@log_func_runtime
def create_expression_summary_cell_type_cube(corpus_path: str):
    """
    Create the cell type expression summary cube. The cell type expression summary cube contains the same
    data as the default expression summary cube, just with a different schema.
    """
    pipeline_state = load_pipeline_state(corpus_path=corpus_path)

    if not pipeline_state.get(EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG):
        raise PipelineStepMissing("expression_summary_default")

    logger.info("Creating the cell type expression summary cube.")
    expression_summary_default_uri = os.path.join(corpus_path, EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME)
    expression_summary_cell_type_uri = os.path.join(corpus_path, EXPRESSION_SUMMARY_CELL_TYPE_CUBE_NAME)

    ctx = create_ctx()
    with tiledb.scope_ctx(ctx):
        dfs = []
        with tiledb.open(expression_summary_default_uri, "r") as cube:
            for row in cube.query(return_incomplete=True).df[:]:
                dfs.append(row)
        expression_summary_default_df = pd.concat(dfs, axis=0)

        create_empty_cube_if_needed(expression_summary_cell_type_uri, expression_summary_schema)
        logger.info(f"Writing cube to {expression_summary_cell_type_uri}")
        tiledb.from_pandas(expression_summary_cell_type_uri, expression_summary_default_df, mode="append")

    pipeline_state[EXPRESSION_SUMMARY_CELL_TYPE_CUBE_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)
