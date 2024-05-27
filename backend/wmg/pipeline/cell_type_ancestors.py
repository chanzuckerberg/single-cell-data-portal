import json
import logging
import os

import tiledb

from backend.common.census_cube.data.snapshot import CELL_COUNTS_CUBE_NAME, CELL_TYPE_ANCESTORS_FILENAME
from backend.common.census_cube.utils import ancestors
from backend.wmg.pipeline.constants import (
    CELL_TYPE_ANCESTORS_CREATED_FLAG,
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
)
from backend.wmg.pipeline.errors import PipelineStepMissing
from backend.wmg.pipeline.utils import load_pipeline_state, log_func_runtime, write_pipeline_state

logger = logging.getLogger(__name__)


@log_func_runtime
def create_cell_type_ancestors(corpus_path: str) -> None:
    """ """
    logger.info("Generating dataset metadata file")
    pipeline_state = load_pipeline_state(corpus_path)
    if not pipeline_state.get(EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG):
        raise PipelineStepMissing("cell_counts")

    with tiledb.open(os.path.join(corpus_path, CELL_COUNTS_CUBE_NAME)) as cell_counts_cube:
        cell_counts_df = cell_counts_cube.df[:]

    valid_cell_types = cell_counts_df["cell_type_ontology_term_id"].unique()
    ancestors_dict = {}
    for cell_type_id in valid_cell_types:
        ancestors_dict[cell_type_id] = sorted(set(ancestors(cell_type_id)).intersection(valid_cell_types))

    logger.info("Writing cell type ancestors file")
    with open(f"{corpus_path}/{CELL_TYPE_ANCESTORS_FILENAME}", "w") as f:
        json.dump(ancestors_dict, f)
    pipeline_state[CELL_TYPE_ANCESTORS_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)
