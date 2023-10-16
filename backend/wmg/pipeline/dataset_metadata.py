import json
import logging

from backend.wmg.data.snapshot import DATASET_METADATA_FILENAME
from backend.wmg.data.utils import get_datasets_from_discover_api, log_func_runtime
from backend.wmg.pipeline.constants import DATASET_METADATA_CREATED_FLAG
from backend.wmg.pipeline.utils import load_pipeline_state, write_pipeline_state

logger = logging.getLogger(__name__)


@log_func_runtime
def create_dataset_metadata(*, corpus_path: str) -> None:
    """
    This function generates a dictionary containing metadata for each dataset.
    The metadata includes the dataset id, label, collection id, and collection label.
    The function fetches the datasets from the Discover API and iterates over them to create the metadata dictionary.
    """
    logger.info("Generating dataset metadata file")
    pipeline_state = load_pipeline_state(corpus_path)
    datasets = get_datasets_from_discover_api()
    dataset_dict = {}
    for dataset in datasets:
        dataset_id = dataset["dataset_id"]
        dataset_dict[dataset_id] = dict(
            id=dataset_id,
            label=dataset["title"],
            collection_id=dataset["collection_id"],
            collection_label=dataset["collection_name"],
        )
    logger.info("Writing dataset metadata file")
    with open(f"{corpus_path}/{DATASET_METADATA_FILENAME}", "w") as f:
        json.dump(dataset_dict, f)
    pipeline_state[DATASET_METADATA_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)
