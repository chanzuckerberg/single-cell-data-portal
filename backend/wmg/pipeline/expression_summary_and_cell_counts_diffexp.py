import logging
import os

import pandas as pd
import tiledb

from backend.wmg.data.schemas.cube_schema_diffexp import (
    cell_counts_logical_dims,
    cell_counts_logical_dims_exclude_dataset_id,
    cell_counts_schema,
    expression_summary_schema,
)
from backend.wmg.data.snapshot import (
    CELL_COUNTS_CUBE_NAME,
    CELL_COUNTS_DIFFEXP_CUBE_NAME,
    EXPRESSION_SUMMARY_CUBE_NAME,
    EXPRESSION_SUMMARY_DIFFEXP_CUBE_NAME,
)
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_DIFFEXP_CUBES_CREATED_FLAG,
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
def create_expression_summary_and_cell_counts_diffexp_cubes(corpus_path: str):
    pipeline_state = load_pipeline_state(corpus_path=corpus_path)

    if not pipeline_state.get(EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG):
        raise PipelineStepMissing("expression_summary")

    logger.info("Creating the differential expression expression summary cubes.")
    expression_summary_uri = os.path.join(corpus_path, EXPRESSION_SUMMARY_CUBE_NAME)
    cell_counts_uri = os.path.join(corpus_path, CELL_COUNTS_CUBE_NAME)

    expression_summary_diffexp_uri = os.path.join(corpus_path, EXPRESSION_SUMMARY_DIFFEXP_CUBE_NAME)
    cell_counts_diffexp_uri = os.path.join(corpus_path, CELL_COUNTS_DIFFEXP_CUBE_NAME)
    ctx = create_ctx()
    with tiledb.scope_ctx(ctx):
        # cell counts
        create_empty_cube_if_needed(cell_counts_diffexp_uri, cell_counts_schema)
        with tiledb.open(cell_counts_uri, "r") as cube:
            cell_counts_df = cube.df[:]

        groups_no_dataset_id = cell_counts_df.groupby(cell_counts_logical_dims_exclude_dataset_id).first().index
        group_ids_indexer = pd.Series(range(len(groups_no_dataset_id)), index=groups_no_dataset_id)

        cell_counts_df = cell_counts_df.groupby(cell_counts_logical_dims).sum(numeric_only=True).reset_index()
        groups = cell_counts_df.set_index(cell_counts_logical_dims_exclude_dataset_id).index
        cell_counts_df["group_id"] = group_ids_indexer[groups].values

        logger.info(f"Writing cell_counts diffexp cube to {cell_counts_diffexp_uri}")
        tiledb.from_pandas(cell_counts_diffexp_uri, cell_counts_df, mode="append")

        # expression summary
        create_empty_cube_if_needed(expression_summary_diffexp_uri, expression_summary_schema)
        logger.info(f"Writing expression_summary diffexp cube to {expression_summary_diffexp_uri}")
        with tiledb.open(expression_summary_uri, "r") as cube:
            for row in cube.query(return_incomplete=True).df[:]:
                row = (
                    row.groupby(cell_counts_logical_dims_exclude_dataset_id + ["gene_ontology_term_id"])
                    .sum(numeric_only=True)
                    .reset_index()
                    .set_index(cell_counts_logical_dims_exclude_dataset_id)
                )
                row["group_id"] = group_ids_indexer[row.index].values
                row = row.reset_index(drop=True)

                tiledb.from_pandas(
                    expression_summary_diffexp_uri,
                    row[_get_columns_from_array_schema(expression_summary_schema)],
                    mode="append",
                )

        # consolidate and vacuum
        for uri in [expression_summary_diffexp_uri, cell_counts_diffexp_uri]:
            tiledb.consolidate(uri, ctx=ctx)
            tiledb.vacuum(uri, ctx=ctx)

    pipeline_state[EXPRESSION_SUMMARY_AND_CELL_COUNTS_DIFFEXP_CUBES_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)


def _get_columns_from_array_schema(array_schema: tiledb.ArraySchema):
    dimension_names = [dim.name for dim in array_schema.domain]
    attribute_names = [attr.name for attr in array_schema]
    return dimension_names + attribute_names
