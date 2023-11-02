import cellxgene_census
import tiledbsoma as soma

from backend.common.feature_flag import FeatureFlagService, FeatureFlagValues
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
        organism_mapping = {
            "homo_sapiens": f"is_primary_data == True and nnz >= {GENE_EXPRESSION_COUNT_MIN_THRESHOLD}",
            "mus_musculus": f"is_primary_data == True and nnz >= {GENE_EXPRESSION_COUNT_MIN_THRESHOLD}",
        }
        value_filter = organism_mapping[organism]
        if FeatureFlagService.is_enabled(FeatureFlagValues.SCHEMA_4):
            # Filter out non-tissue tissues and system-level tissues
            # TODO: Once cellxgene-census is updated to support tiledbsoma 1.5.0, we can update the tissue filter to use `not in`:
            # tissue_ontology_term_id not in ['UBERON_0001017', 'UBERON:0001007', 'UBERON:0002405', 'UBERON:0000990', 'UBERON:0001004', 'UBERON:0001434']
            # This will require updating the pinned versions of cellxgene-census and tiledbsoma in `requirements-wmg-pipeline.txt`
            value_filter += " and tissue_type == 'tissue' and tissue_general_ontology_term_id != 'UBERON:0001017' and tissue_general_ontology_term_id != 'UBERON:0001007' and tissue_general_ontology_term_id != 'UBERON:0002405' and tissue_general_ontology_term_id != 'UBERON:0000990' and tissue_general_ontology_term_id != 'UBERON:0001004' and tissue_general_ontology_term_id != 'UBERON:0001434'"
        return value_filter


def create_expression_summary_and_cell_counts_cubes(corpus_path: str):
    pipeline_state = load_pipeline_state(corpus_path)

    with cellxgene_census.open_soma(
        census_version=CensusParameters.census_version,
    ) as census:
        for organismInfo in ORGANISM_INFO:
            organism = organismInfo["label"]
            organismId = organismInfo["id"]

            value_filter = CensusParameters.value_filter(organism)
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
