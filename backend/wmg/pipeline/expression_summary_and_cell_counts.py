import os

import cellxgene_census
import tiledb
import tiledbsoma as soma
from packaging import version
from packaging.specifiers import SpecifierSet

from backend.common.census_cube.data.snapshot import (
    CELL_COUNTS_CUBE_NAME,
)
from backend.wmg.pipeline.cell_counts import create_cell_counts_cube
from backend.wmg.pipeline.constants import (
    ADMISSIBLE_CENSUS_SCHEMA_VERSION_SPECIFIER_SET,
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    ORGANISM_INFO,
    CensusParameters,
)
from backend.wmg.pipeline.expression_summary import ExpressionSummaryCubeBuilder
from backend.wmg.pipeline.utils import (
    load_pipeline_state,
    write_pipeline_state,
)


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

    with cellxgene_census.open_soma(census_version=CensusParameters.census_version) as census:
        census_schema_version, census_build_date = get_census_version_and_build_date(census)

        # Validate the census schema version against an admissible SpecifierSet.
        specifier_set = SpecifierSet(ADMISSIBLE_CENSUS_SCHEMA_VERSION_SPECIFIER_SET)
        parsed_version = version.parse(census_schema_version)
        if parsed_version not in specifier_set:
            raise ValueError(
                f"Unsupported CELLxGENE Census schema version: {census_schema_version}. "
                f"Please use a version of CELLxGENE Census with schema version {ADMISSIBLE_CENSUS_SCHEMA_VERSION_SPECIFIER_SET}."
            )

        dataset_metadata = census["census_info"]["datasets"].read().concat().to_pandas()

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
                        dataset_metadata=dataset_metadata, query=query, corpus_path=corpus_path, organismId=organismId
                    ).create_expression_summary_cube()
                    create_cell_counts_cube(
                        dataset_metadata=dataset_metadata, query=query, corpus_path=corpus_path, organismId=organismId
                    )

    # write census schema version to cell counts cube metadata
    cell_counts_uri = os.path.join(corpus_path, CELL_COUNTS_CUBE_NAME)
    with tiledb.open(cell_counts_uri, mode="w") as A:
        A.meta["census_schema_version"] = census_schema_version
        A.meta["census_build_date"] = census_build_date

    pipeline_state[EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)
