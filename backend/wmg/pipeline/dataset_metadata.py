import json
import logging
import os

import cellxgene_census
import tiledb

from backend.common.census_cube.data.snapshot import CELL_COUNTS_CUBE_NAME, DATASET_METADATA_FILENAME
from backend.wmg.pipeline.constants import (
    DATASET_METADATA_CREATED_FLAG,
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    CensusParameters,
)
from backend.wmg.pipeline.errors import PipelineStepMissing
from backend.wmg.pipeline.utils import load_pipeline_state, log_func_runtime, write_pipeline_state

logger = logging.getLogger(__name__)


@log_func_runtime
def create_dataset_metadata(corpus_path: str) -> None:
    """
    This function generates a dictionary containing metadata for each dataset.
    The metadata includes the dataset id, label, collection id, and collection label.
    The function fetches the datasets from the Discover API and iterates over them to create the metadata dictionary.
    """
    logger.info("Generating dataset metadata file")
    pipeline_state = load_pipeline_state(corpus_path)

    if not pipeline_state.get(EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG):
        raise PipelineStepMissing("cell_counts")

    with cellxgene_census.open_soma(census_version=CensusParameters.census_version) as census:
        dataset_metadata = census["census_info"]["datasets"].read().concat().to_pandas()

    # read in the cell counts df and only keep the dataset_ids that are in the cube
    with tiledb.open(os.path.join(corpus_path, CELL_COUNTS_CUBE_NAME)) as cc_cube:
        cell_counts_df = cc_cube.df[:]
        unique_dataset_ids = cell_counts_df["dataset_id"].unique()

    dataset_metadata = dataset_metadata[dataset_metadata["dataset_id"].isin(unique_dataset_ids)]

    datasets = dataset_metadata.to_dict(orient="records")

    dataset_dict = {}
    for dataset in datasets:
        dataset_dict[dataset["dataset_id"]] = dict(
            id=dataset["dataset_id"],
            label=dataset["dataset_title"],
            collection_id=dataset["collection_id"],
            collection_label=dataset["collection_name"],
        )

    logger.info("Writing dataset metadata file")
    with open(f"{corpus_path}/{DATASET_METADATA_FILENAME}", "w") as f:
        json.dump(dataset_dict, f)
    pipeline_state[DATASET_METADATA_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)
