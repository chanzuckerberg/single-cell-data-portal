import unicodedata

import numba as nb
import numpy as np
import pandas as pd

from backend.wmg.data.utils import (
    get_collections_from_curation_api,
    get_datasets_from_curation_api,
    to_dict,
)


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
    Xs = []
    Ys = []
    for i in range(mat.shape[0]):
        Xs.extend(np.repeat(mat[i], mat[i].size))
        Ys.extend(np.tile(mat[i], mat[i].size))

    # get all the unique combinations of filters
    Xs, Ys = np.unique(np.array((Xs, Ys)), axis=1)

    # exclude self-relationships
    filt = Xs != Ys
    Xs, Ys = Xs[filt], Ys[filt]

    # convert the edges to a linked list representation
    filter_relationships_linked_list = to_dict(Xs, Ys)

    # reorganize the linked list representation to a nested linked list representation
    # where the filter columns are separated into distinct dictionaries
    # e.g. instead of {"cell_type_ontology_term_id__beta cell": ["dataset_id__Single cell transcriptome analysis of human pancreas", "assay_ontology_term_id__assay_type", ...], ...},
    # it's now {"cell_type_ontology_term_id__beta cell": {"dataset_id": ["dataset_id__Single cell transcriptome analysis of human pancreas", ...], "assay_ontology_term_id": ["assay_ontology_term_id__assay_type", ...], ...}, ...}.
    # This structure is easier to parse by the `/query` endpoint.
    for k, v in filter_relationships_linked_list.items():
        filter_relationships_linked_list[k] = to_dict([x.split("__")[0] for x in v], v)

    return filter_relationships_linked_list


@nb.njit(fastmath=True, error_model="numpy", parallel=False, nogil=True)
def gene_expression_sum_x_cube_dimension(
    rankit_values: np.ndarray,
    obs_idxs: np.ndarray,
    var_idx: np.ndarray,
    cube_indices: np.ndarray,
    sum_into: np.ndarray,
    sqsum_into: np.ndarray,
    nnz_into: np.ndarray,
):
    """
    Sum the rankit values for each gene (for each cube row/combo of cell attributes)
    Also track the number of cells that express that gene (nnz count)
    """
    for k in range(len(rankit_values)):
        val = rankit_values[k]
        if np.isfinite(val):
            cidx = var_idx[k]
            grp_idx = cube_indices[obs_idxs[k]]
            sum_into[grp_idx, cidx] += val
            sqsum_into[grp_idx, cidx] += val**2
            nnz_into[grp_idx, cidx] += 1


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
