import logging
import os

import pandas as pd
import tiledb

from backend.wmg.data.schemas.expression_summary_cube_schemas_diffexp import (
    base_expression_summary_indexed_dims,
    expression_summary_attrs,
    expression_summary_schemas,
)
from backend.wmg.data.snapshot import (
    EXPRESSION_SUMMARY_CUBE_NAME,
    EXPRESSION_SUMMARY_DIFFEXP_CUBE_PREFIX,
)
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    EXPRESSION_SUMMARY_DIFFEXP_CUBES_CREATED_FLAG,
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
def create_expression_summary_diffexp_cubes(corpus_path: str):
    pipeline_state = load_pipeline_state(corpus_path=corpus_path)

    if not pipeline_state.get(EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG):
        raise PipelineStepMissing("expression_summary")

    logger.info("Creating the default expression summary cube.")
    expression_summary_uri = os.path.join(corpus_path, EXPRESSION_SUMMARY_CUBE_NAME)

    cube_uris = {
        secondary_dim: os.path.join(corpus_path, f"{EXPRESSION_SUMMARY_DIFFEXP_CUBE_PREFIX}_{secondary_dim}")
        for secondary_dim in expression_summary_schemas
    }

    ctx = create_ctx()
    with tiledb.scope_ctx(ctx):
        dfs = {secondary_dim: [] for secondary_dim in cube_uris}
        with tiledb.open(expression_summary_uri, "r") as cube:
            for row in cube.query(return_incomplete=True).df[:]:
                # we are generating aggregated dataframes online to avoid reading through the original cube
                # multiple times. This is a tradeoff between memory and time.
                # Note: if memory becomes an issue, we can switch to performing the aggregation one secondary
                # dimension at a time.
                for secondary_dim in cube_uris:
                    groupby_dims = base_expression_summary_indexed_dims + expression_summary_attrs + [secondary_dim]
                    dfs[secondary_dim].append(row.groupby(groupby_dims).sum(numeric_only=True).reset_index())

        for secondary_dim, uri in cube_uris.items():
            create_empty_cube_if_needed(uri, expression_summary_schemas[secondary_dim])
            groupby_dims = base_expression_summary_indexed_dims + expression_summary_attrs + [secondary_dim]
            expression_summary_df = (
                pd.concat(dfs[secondary_dim], axis=0).groupby(groupby_dims).sum(numeric_only=True).reset_index()
            )

            logger.info(f"Writing cube to {uri}")
            tiledb.from_pandas(uri, pd.concat(expression_summary_df, axis=0), mode="append")

    pipeline_state[EXPRESSION_SUMMARY_DIFFEXP_CUBES_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)
