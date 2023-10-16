import logging
import pathlib
import time

from backend.wmg.data.utils import log_func_runtime
from backend.wmg.pipeline.cell_type_ordering import create_cell_type_ordering
from backend.wmg.pipeline.constants import (
    CELL_TYPE_ORDERING_CREATED_FLAG,
    DATASET_METADATA_CREATED_FLAG,
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG,
    FILTER_RELATIONSHIPS_CREATED_FLAG,
    MARKER_GENES_CUBE_CREATED_FLAG,
    PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG,
)
from backend.wmg.pipeline.dataset_metadata import create_dataset_metadata
from backend.wmg.pipeline.expression_summary_and_cell_counts import create_expression_summary_and_cell_counts_cubes
from backend.wmg.pipeline.expression_summary_default import create_expression_summary_default_cube
from backend.wmg.pipeline.filter_relationships import create_filter_relationships_graph
from backend.wmg.pipeline.marker_genes import create_marker_genes_cube
from backend.wmg.pipeline.primary_filter_dimensions import create_primary_filter_dimensions
from backend.wmg.pipeline.utils import load_pipeline_state

logger = logging.getLogger(__name__)

PIPELINE_STEPS = [
    {
        "flag": EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
        "step": create_expression_summary_and_cell_counts_cubes,
    },
    {"flag": EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG, "step": create_expression_summary_default_cube},
    {"flag": FILTER_RELATIONSHIPS_CREATED_FLAG, "step": create_filter_relationships_graph},
    {"flag": PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG, "step": create_primary_filter_dimensions},
    {"flag": MARKER_GENES_CUBE_CREATED_FLAG, "step": create_marker_genes_cube},
    {"flag": DATASET_METADATA_CREATED_FLAG, "step": create_dataset_metadata},
    {"flag": CELL_TYPE_ORDERING_CREATED_FLAG, "step": create_cell_type_ordering},
]


@log_func_runtime
def run_pipeline():
    corpus_path = str(int(time.time()))
    logger.info(f"Creating directory {corpus_path}")
    pathlib.Path(corpus_path).mkdir(parents=True, exist_ok=True)

    pipeline_state = load_pipeline_state(corpus_path)

    for pipeline_step in PIPELINE_STEPS:
        flag = pipeline_step["flag"]
        step = pipeline_step["step"]

        if not pipeline_state.get(flag):
            step(corpus_path=corpus_path)
