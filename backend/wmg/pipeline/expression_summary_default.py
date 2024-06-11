import logging
import os

import pandas as pd
import tiledb

from backend.common.census_cube.data.schemas.cube_schema_default import (
    expression_summary_indexed_dims,
    expression_summary_non_indexed_dims,
    expression_summary_schema,
)
from backend.common.census_cube.data.snapshot import (
    EXPRESSION_SUMMARY_CUBE_NAME,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME,
)
from backend.common.census_cube.data.tiledb import create_ctx
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
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
def create_expression_summary_default_cube(corpus_path: str):
    """
    Create the default expression summary cube. The default expression summary cube is an aggregation across
    non-default dimensions in the expression summary cube.
    """
    pipeline_state = load_pipeline_state(corpus_path=corpus_path)

    if not pipeline_state.get(EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG):
        raise PipelineStepMissing("expression_summary")

    logger.info("Creating the default expression summary cube.")
    expression_summary_uri = os.path.join(corpus_path, EXPRESSION_SUMMARY_CUBE_NAME)
    expression_summary_default_uri = os.path.join(corpus_path, EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME)

    ctx = create_ctx()
    with tiledb.scope_ctx(ctx):
        dfs = []
        with tiledb.open(expression_summary_uri, "r") as cube:
            for row in cube.query(return_incomplete=True).df[:]:
                dfs.append(row)
        expression_summary_df = pd.concat(dfs, axis=0)

        # TODO: Explore if staying in PyArrow space incurs less memory overhead. If so, use PyArrow groupby
        expression_summary_df_default = (
            expression_summary_df.groupby(expression_summary_indexed_dims + expression_summary_non_indexed_dims)
            .sum(numeric_only=True)
            .reset_index()
        )

        create_empty_cube_if_needed(expression_summary_default_uri, expression_summary_schema)
        logger.info(f"Writing cube to {expression_summary_default_uri}")
        tiledb.from_pandas(expression_summary_default_uri, expression_summary_df_default, mode="append")

    pipeline_state[EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)
