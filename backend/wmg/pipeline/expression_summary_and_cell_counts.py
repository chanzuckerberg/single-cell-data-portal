import cellxgene_census
import tiledbsoma as soma

from backend.wmg.pipeline.cell_counts import create_cell_counts_cube
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    GENE_EXPRESSION_COUNT_MIN_THRESHOLD,
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


class CensusParameters:
    census_version = "latest"

    def value_filter(organism: str) -> str:
        value_filter_mapping = {
            "homo_sapiens": f"is_primary_data == True and nnz >= {GENE_EXPRESSION_COUNT_MIN_THRESHOLD}",
            "mus_musculus": f"is_primary_data == True and nnz >= {GENE_EXPRESSION_COUNT_MIN_THRESHOLD}",
        }
        value_filter = value_filter_mapping[organism]
        # if schema_4:
        #     value_filter += " and tissue_type == 'tissue'"
        return value_filter


def create_expression_summary_and_cell_counts_cubes(corpus_path: str):
    pipeline_state = load_pipeline_state(corpus_path)

    with cellxgene_census.open_soma(
        census_version=CensusParameters.census_version,
    ) as census:
        for organismInfo in ORGANISM_INFO:
            organism = organismInfo["label"]
            organismId = organismInfo["id"]

            value_filter = CensusParameters.value_filter[organism]
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
