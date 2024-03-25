import os

import cellxgene_census
import tiledb
import tiledbsoma as soma
from packaging import version

from backend.wmg.data.snapshot import (
    CELL_COUNTS_CUBE_NAME,
)
from backend.wmg.pipeline.cell_counts import create_cell_counts_cube
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    GENE_EXPRESSION_COUNT_MIN_THRESHOLD,
    MAXIMUM_ADMISSIBLE_CENSUS_SCHEMA_MAJOR_VERSION,
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
        # Filter out system-level tissues. Census filters out organoids + cell cultures
        value_filter += " and tissue_general_ontology_term_id != 'UBERON:0001017' and tissue_general_ontology_term_id != 'UBERON:0001007' and tissue_general_ontology_term_id != 'UBERON:0002405' and tissue_general_ontology_term_id != 'UBERON:0000990' and tissue_general_ontology_term_id != 'UBERON:0001004' and tissue_general_ontology_term_id != 'UBERON:0001434'"
        return value_filter


def get_census_version_and_build_date(census: soma.Collection):
    """
    Retrieves the census schema version and build date from the given census collection.

    Parameters:
    census (soma.Collection): The census collection to retrieve the schema version and build date from.

    Returns:
    tuple: A tuple containing the census schema version and build date.
    """
    census_schema_version, census_build_date = (
        census["census_info"]["summary"]
        .read()
        .concat()
        .to_pandas()
        .set_index("label")
        .loc[["census_schema_version", "census_build_date"], "value"]
    ).tolist()

    return census_schema_version, census_build_date


def create_expression_summary_and_cell_counts_cubes(corpus_path: str):
    pipeline_state = load_pipeline_state(corpus_path)
    context = soma.options.SOMATileDBContext()
    tiledb_config = {}
    tiledb_config["vfs.s3.region"] = "us-west-2"
    # S3 requests should not be signed, since we want to allow anonymous access
    tiledb_config["vfs.s3.no_sign_request"] = "false"
    context = context.replace(tiledb_config=tiledb_config)
    with cellxgene_census.open_soma(
        uri="s3://bruce-tmp/census-schema-five-test-build-2/soma/", context=context
    ) as census:
        census_schema_version, census_build_date = get_census_version_and_build_date(census)

        major_census_schema_version = version.parse(census_schema_version).major
        if major_census_schema_version > MAXIMUM_ADMISSIBLE_CENSUS_SCHEMA_MAJOR_VERSION:
            raise ValueError(
                f"Unsupported census schema version: {census_schema_version}. "
                f"Please use a version of cellxgene-census that supports census schema version {MAXIMUM_ADMISSIBLE_CENSUS_SCHEMA_MAJOR_VERSION} or lower."
            )

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

    # write census schema version to cell counts cube metadata
    cell_counts_uri = os.path.join(corpus_path, CELL_COUNTS_CUBE_NAME)
    with tiledb.open(cell_counts_uri, mode="w") as A:
        A.meta["census_schema_version"] = census_schema_version
        A.meta["census_build_date"] = census_build_date

    pipeline_state[EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)
