import json
import logging
import os
import unicodedata

from tiledb import ArraySchema

from backend.wmg.data.tissue_mapper import TissueMapper
from backend.wmg.data.utils import (
    create_empty_cube,
    get_collections_from_curation_api,
    get_datasets_from_curation_api,
)
from backend.wmg.pipeline.constants import (
    DATASET_METADATA_CREATED_FLAG,
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG,
    FILTER_RELATIONSHIPS_CREATED_FLAG,
    MARKER_GENES_CUBE_CREATED_FLAG,
    PIPELINE_STATE_FILENAME,
    PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG,
)

logger = logging.getLogger(__name__)


def load_pipeline_state(corpus_path: str):
    state_file = os.path.join(corpus_path, PIPELINE_STATE_FILENAME)
    if os.path.exists(state_file):
        with open(state_file, "r") as f:
            state = json.load(f)
    else:
        state = {
            EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG: False,
            EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG: False,
            MARKER_GENES_CUBE_CREATED_FLAG: False,
            FILTER_RELATIONSHIPS_CREATED_FLAG: False,
            PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG: False,
            DATASET_METADATA_CREATED_FLAG: False,
        }
    return state


def write_pipeline_state(pipeline_state: dict, corpus_path: str):
    with open(os.path.join(corpus_path, PIPELINE_STATE_FILENAME), "w") as f:
        json.dump(pipeline_state, f)


def remove_accents(input_str):
    nfkd_form = unicodedata.normalize("NFKD", input_str)
    return "".join([c for c in nfkd_form if not unicodedata.combining(c)])


def return_dataset_dict_w_publications():
    datasets = get_datasets_from_curation_api()
    collections = get_collections_from_curation_api()
    collections_dict = {collection["collection_id"]: collection for collection in collections}

    # cchoi: creating a helper function to format citations properly
    def create_formatted_citation(collection):
        publisher_metadata = collection["publisher_metadata"]
        if publisher_metadata is None:
            return "No Publication"
        first_author = collection["publisher_metadata"]["authors"][0]
        # first_author could be either 'family' or 'name'
        citation = f"{first_author['family'] if 'family' in first_author else first_author['name']} et al. {collection['publisher_metadata']['journal']} {collection['publisher_metadata']['published_year']}"
        formatted_citation = "No Publication" if collection["publisher_metadata"]["is_preprint"] else citation
        return formatted_citation

    dataset_dict = {}
    for dataset in datasets:
        dataset_id = dataset["dataset_id"]
        collection = collections_dict[dataset["collection_id"]]
        dataset_dict[dataset_id] = create_formatted_citation(collection)

    return dataset_dict


def list_grouped_primary_filter_dimensions_term_ids(
    df, primary_dim_name: str, group_by_dim: str
) -> dict[str, list[str]]:
    """
    This function takes a dataframe and two dimension names as input. It groups the dataframe by the second dimension name,
    and lists the unique values of the first dimension name for each group. The output is a dictionary where the keys are the
    unique values of the second dimension name, and the values are lists of unique values of the first dimension name for each group.

    :param df: The input dataframe.
    :param primary_dim_name: The name of the first dimension.
    :param group_by_dim: The name of the second dimension to group by.
    :return: A dictionary of lists of unique values of the first dimension for each group.
    """
    return (
        df[[primary_dim_name, group_by_dim]]
        .drop_duplicates()
        .groupby(group_by_dim)
        .agg(list)
        .to_dict()[primary_dim_name]
    )


def order_tissues(ontology_term_ids: list[str]) -> list[str]:
    """
    Order tissues based on appearance in TissueMapper.HIGH_LEVEL_TISSUES. This will maintain the priority set in
    that class which is intended to keep most relevant tissues on top and tissues that are related to be placed
    sequentially
    """
    ontology_term_ids = set(ontology_term_ids)
    ordered_ontology_term_ids = []
    for tissue in TissueMapper.HIGH_LEVEL_TISSUES:
        tissue = TissueMapper.reformat_ontology_term_id(tissue, to_writable=True)
        if tissue in ontology_term_ids:
            ontology_term_ids.remove(tissue)
            ordered_ontology_term_ids.append(tissue)

    if ontology_term_ids:
        ordered_ontology_term_ids = ordered_ontology_term_ids + list(ontology_term_ids)

    return ordered_ontology_term_ids


def create_empty_cube_if_needed(uri: str, schema: ArraySchema):
    if not os.path.exists(uri):
        logger.info(f"Creating empty cube at {uri} with schema: {schema}")
        create_empty_cube(uri, schema)
