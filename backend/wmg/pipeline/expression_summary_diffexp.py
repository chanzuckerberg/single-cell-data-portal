import json
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
    CARDINALITY_PER_DIMENSION_FILENAME,
    EXPRESSION_SUMMARY_CUBE_NAME,
    EXPRESSION_SUMMARY_DIFFEXP_CUBE_NAMES,
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
        cube_name.split("__")[-1]: os.path.join(corpus_path, cube_name)
        for cube_name in EXPRESSION_SUMMARY_DIFFEXP_CUBE_NAMES
    }

    default_groupby_dims = base_expression_summary_indexed_dims + expression_summary_attrs
    ctx = create_ctx()
    with tiledb.scope_ctx(ctx):

        for secondary_dim, uri in cube_uris.items():
            create_empty_cube_if_needed(uri, expression_summary_schemas[secondary_dim])

        expression_summary_default_df_chunks = []
        cardinality_per_dimension = {}
        with tiledb.open(expression_summary_uri, "r") as cube:
            for row in cube.query(return_incomplete=True).df[:]:
                # write the row chunk into the new cube with its specific ArraySchema
                for secondary_dim, uri in cube_uris.items():
                    if secondary_dim == "default":
                        continue

                    # update the cardinality per dimension
                    row_cats = row.select_dtypes(exclude="number")
                    for col in row_cats.columns:
                        L = cardinality_per_dimension.get(col, [])
                        L.extend(row_cats[col].unique())
                        L = list(set(L))
                        cardinality_per_dimension[col] = L

                    tiledb.from_pandas(
                        uri,
                        row[_get_columns_from_array_schema(expression_summary_schemas[secondary_dim])],
                        mode="append",
                    )

                # generate the default cube chunk by chunk
                expression_summary_default_df_chunks.append(
                    row.groupby(default_groupby_dims).sum(numeric_only=True).reset_index()
                )

            # perform a final groupby and write the default cube
            expression_summary_default_df = (
                pd.concat(expression_summary_default_df_chunks, axis=0)
                .groupby(default_groupby_dims)
                .sum(numeric_only=True)
                .reset_index()
            )

            cardinality_per_dimension = {col: len(L) for col, L in cardinality_per_dimension.items()}

            logger.info(f"Writing cube to {cube_uris['default']}")
            tiledb.from_pandas(
                cube_uris["default"],
                expression_summary_default_df[_get_columns_from_array_schema(expression_summary_schemas["default"])],
                mode="append",
            )

            with open(os.path.join(corpus_path, CARDINALITY_PER_DIMENSION_FILENAME), "w") as f:
                json.dump(cardinality_per_dimension, f)

    pipeline_state[EXPRESSION_SUMMARY_DIFFEXP_CUBES_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)


def _get_columns_from_array_schema(array_schema: tiledb.ArraySchema):
    dimension_names = [dim.name for dim in array_schema.domain]
    attribute_names = [attr.name for attr in array_schema]
    return dimension_names + attribute_names
