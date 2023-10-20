import os

import cellxgene_census
import tiledbsoma as soma

from backend.wmg.pipeline.cell_counts import create_cell_counts_cube
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    GENE_EXPRESSION_COUNT_MIN_THRESHOLD,
    WMG_PIPELINE_TEST_RUN_CONFIG,
    WMG_PIPELINE_TEST_RUN_KEY,
)
from backend.wmg.pipeline.expression_summary import ExpressionSummaryCubeBuilder
from backend.wmg.pipeline.utils import (
    load_pipeline_state,
    write_pipeline_state,
)

ORGANISM_INFO = [
    {"label": "homo_sapiens", "id": "NCBITaxon:9606"},
    {"label": "mus_musculus", "id": "NCBITaxon:10090"},
]


def create_expression_summary_and_cell_counts_cubes(corpus_path: str):
    pipeline_state = load_pipeline_state(corpus_path)

    WMG_PIPELINE_TEST_RUN = os.environ.get(WMG_PIPELINE_TEST_RUN_KEY) == "true"
    for organismInfo in ORGANISM_INFO:
        organism = organismInfo["label"]
        organismId = organismInfo["id"]

        value_filter = (
            WMG_PIPELINE_TEST_RUN_CONFIG["value_filters"][organism]
            if WMG_PIPELINE_TEST_RUN
            else f"is_primary_data == True and nnz >= {GENE_EXPRESSION_COUNT_MIN_THRESHOLD}"
        )
        census_version = WMG_PIPELINE_TEST_RUN_CONFIG["census_version"] if WMG_PIPELINE_TEST_RUN else "latest"

        with cellxgene_census.open_soma(
            census_version=census_version,
        ) as census:
            organism = census["census_data"][organism]
            with organism.axis_query(
                "RNA",
                obs_query=soma.AxisQuery(value_filter=value_filter),
            ) as query:
                if not pipeline_state.get(EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG):
                    ExpressionSummaryCubeBuilder(
                        query=query, corpus_path=corpus_path, organismId=organismId
                    ).create_expression_summary_cube()
                    create_cell_counts_cube(query=query, corpus_path=corpus_path, organismId=organismId)

    pipeline_state[EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)
