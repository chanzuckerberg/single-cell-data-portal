import json
import logging
import unicodedata

import numpy as np
import pandas as pd
import tiledb

from backend.wmg.data.schemas.corpus_schema import (
    FILTER_RELATIONSHIPS_NAME,
    OBS_ARRAY_NAME,
)
from backend.wmg.data.schemas.cube_schema import cell_counts_schema
from backend.wmg.data.snapshot import CELL_COUNTS_CUBE_NAME
from backend.wmg.data.utils import (
    create_empty_cube,
    get_collections_from_curation_api,
    get_datasets_from_curation_api,
    log_func_runtime,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def extract(corpus_path: str) -> pd.DataFrame:
    """
    get obs data from integrated corpus
    """
    tiledb_array = tiledb.open(f"{corpus_path}/{OBS_ARRAY_NAME}")
    return tiledb_array.df[:]


def transform(obs: pd.DataFrame) -> pd.DataFrame:
    """
    Create cell count cube data frame by grouping data in the
    integrated corpus obs arrays on relevant features
    """

    # filter out observations in the 'filter_cells' attribute.
    # It is important to filter out rejected observations BEFORE
    # performing the `groupby` operation because the `groupby` operation
    # would lose information about the `filter_cells` attribute because
    # `filter_cells` is not one of the columns used in the `groupby` list
    obs_to_keep = obs[np.logical_not(obs["filter_cells"])]

    df = (
        obs_to_keep.groupby(
            by=[
                "dataset_id",
                "cell_type_ontology_term_id",
                "tissue_ontology_term_id",
                "tissue_original_ontology_term_id",
                "assay_ontology_term_id",
                "development_stage_ontology_term_id",
                "disease_ontology_term_id",
                "self_reported_ethnicity_ontology_term_id",
                "sex_ontology_term_id",
                "organism_ontology_term_id",
            ],
            as_index=False,
        ).size()
    ).rename(columns={"size": "n_cells"})

    return df


def load(corpus_path: str, df: pd.DataFrame) -> str:
    """
    write cell count cube to disk
    """
    uri = f"{corpus_path}/{CELL_COUNTS_CUBE_NAME}"
    create_empty_cube(uri, cell_counts_schema)
    tiledb.from_pandas(uri, df, mode="append")
    return uri


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


@log_func_runtime
def create_cell_count_cube(corpus_path: str):
    """
    Create cell count cube and write to disk
    """
    obs = extract(corpus_path)
    df = transform(obs)

    dataset_dict = return_dataset_dict_w_publications()
    df["publication_citation"] = [
        remove_accents(dataset_dict.get(dataset_id, "No Publication")) for dataset_id in df["dataset_id"]
    ]

    n_cells = df["n_cells"].to_numpy()
    df["n_cells"] = n_cells

    filter_relationships_linked_list = create_filter_relationships_graph(df)

    with open(f"{corpus_path}/{FILTER_RELATIONSHIPS_NAME}.json", "w") as f:
        json.dump(filter_relationships_linked_list, f)

    uri = load(corpus_path, df)
    cell_count = df.n_cells.sum()
    logger.info(f"{cell_count=}")
    logger.info(f"Cell count cube created and stored at {uri}")


def create_filter_relationships_graph(df: pd.DataFrame) -> dict:
    """
    Create a graph of filter relationships

    Arguments
    ---------
    df: pd.DataFrame
        Dataframe containing the filter combinations per cell

    Returns
    -------
    filter_relationships_hash_table: dict
        A dictionary containing the filter relationships for each unique filter value.
        Filter values are their column name + __ + their value (e.g. "dataset_id__Single cell transcriptome analysis of human pancreas").
        The dictionary values are lists of filter values that are related to the key filter value. Relatedness indicates that these filters
        are co-occuring in at least one cell.
    """
    # get a dataframe of the columns that are not numeric
    df_filters = df.select_dtypes(exclude="number")
    # get a numpy array of the column names with shape (1, n_cols)
    cols = df_filters.columns.values[None, :]

    # tile the column names row-wise to match the shape of the dataframe and concatenate to the values
    # this ensures that filter values will never collide across columns.
    mat = np.tile(cols, (df.shape[0], 1)) + "__" + df_filters.values

    # for each cell, get all pairwise combinations of filters compresent in that cell
    # these are the edges of the filter relationships graph
    Xs = []  # type: ignore
    Ys = []  # type: ignore
    for i in range(mat.shape[0]):
        Xs.extend(np.repeat(mat[i], mat[i].size))
        Ys.extend(np.tile(mat[i], mat[i].size))

    # get all the unique combinations of filters
    Xs, Ys = np.unique(np.array((Xs, Ys)), axis=1)

    # exclude self-relationships
    filt = Xs != Ys
    Xs, Ys = Xs[filt], Ys[filt]

    # convert the edges to a linked list representation
    filter_relationships_linked_list = _to_dict(Xs, Ys)

    # reorganize the linked list representation to a nested linked list representation
    # where the filter columns are separated into distinct dictionaries
    # e.g. instead of {"cell_type_ontology_term_id__beta cell": ["dataset_id__Single cell transcriptome analysis of human pancreas", "assay_ontology_term_id__assay_type", ...], ...},
    # it's now {"cell_type_ontology_term_id__beta cell": {"dataset_id": ["dataset_id__Single cell transcriptome analysis of human pancreas", ...], "assay_ontology_term_id": ["assay_ontology_term_id__assay_type", ...], ...}, ...}.
    # This structure is easier to parse by the `/query` endpoint.
    for k, v in filter_relationships_linked_list.items():
        filter_relationships_linked_list[k] = _to_dict([x.split("__")[0] for x in v], v)

    return filter_relationships_linked_list


def _to_dict(a, b):
    """
    convert a flat key array (a) and a value array (b) into a dictionary with values grouped by keys
    """
    a = np.array(a)
    b = np.array(b)
    idx = np.argsort(a)
    a = a[idx]
    b = b[idx]
    bounds = np.where(a[:-1] != a[1:])[0] + 1
    bounds = np.append(np.append(0, bounds), a.size)
    bounds_left = bounds[:-1]
    bounds_right = bounds[1:]
    slists = [b[bounds_left[i] : bounds_right[i]] for i in range(bounds_left.size)]
    d = dict(zip(np.unique(a), [list(set(x)) for x in slists]))
    return d


def remove_accents(input_str):
    nfkd_form = unicodedata.normalize("NFKD", input_str)
    return "".join([c for c in nfkd_form if not unicodedata.combining(c)])
