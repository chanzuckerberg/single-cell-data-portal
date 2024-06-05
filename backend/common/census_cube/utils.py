import os
from functools import lru_cache
from typing import Dict, Optional

import numba as nb
import numpy as np
import pandas as pd
import requests
from cellxgene_ontology_guide.ontology_parser import OntologyParser
from requests.adapters import HTTPAdapter
from urllib3.util import Retry

from backend.common.census_cube.data.constants import CENSUS_CUBE_PINNED_SCHEMA_VERSION
from backend.common.census_cube.data.snapshot import CensusCubeSnapshot
from backend.common.constants import DEPLOYMENT_STAGE_TO_API_URL

# exported and used by all modules related to the census cube
ontology_parser = OntologyParser(schema_version=f"v{CENSUS_CUBE_PINNED_SCHEMA_VERSION}")


def find_all_dim_option_values(snapshot, organism: str, dimension: str) -> list:
    organism_key = "organism_ontology_term_id__" + organism
    all_filter_options = snapshot.filter_relationships[organism_key].get(dimension, [])
    return [option.split("__")[1] for option in all_filter_options]


def find_dim_option_values(criteria: Dict, snapshot, dimension: str) -> list:
    """Find values for the specified dimension that satisfy the given filtering criteria,
    ignoring any criteria specified for the given dimension."""

    filter_options_criteria = dict(criteria)
    # Remove gene_ontology_term_ids from the criteria as it is not an eligible cross-filter dimension.
    filter_options_criteria.pop("gene_ontology_term_ids", None)

    # depluralize `dimension` if necessary
    dimension = depluralize(dimension)

    # each element  in `linked_filter_sets` corresponds to the set of filters linked to the attributes specified for a corresponding criteria key
    linked_filter_sets = []

    # `all_criteria_attributes` is the set of all attributes specified across all criteria
    all_criteria_attributes = set()

    for key in filter_options_criteria:
        attrs = filter_options_criteria[key]

        # depluralize `key` if necessary
        key = depluralize(key)

        # ignore the criteria for the specified dimension
        if key != dimension:
            if isinstance(attrs, list):
                if len(attrs) > 0:
                    # prepend the key to each attribute value
                    prefixed_attributes = [key + "__" + val for val in attrs]
                    all_criteria_attributes = all_criteria_attributes.union(prefixed_attributes)

                    # for each attribute (attr) in `prefixed_attributes`,
                    # get the set of filters for the specified dimension that are linked to `attr`
                    linked_filter_set = set()
                    for attr in prefixed_attributes:
                        if dimension in snapshot.filter_relationships.get(attr, {}):
                            linked_filter_set = linked_filter_set.union(
                                set(snapshot.filter_relationships[attr][dimension])
                            )

                    linked_filter_sets.append(linked_filter_set)
            else:
                if attrs != "":
                    prefixed_attribute = key + "__" + attrs
                    all_criteria_attributes.add(prefixed_attribute)
                    if dimension in snapshot.filter_relationships.get(prefixed_attribute, {}):
                        linked_filter_sets.append(set(snapshot.filter_relationships[prefixed_attribute][dimension]))

    # the candidate options are the intersection of the sets of linked filters for each criteria key
    if len(linked_filter_sets) > 1:
        candidate_options = linked_filter_sets[0].intersection(*linked_filter_sets[1:])
    else:
        candidate_options = linked_filter_sets[0]

    # each valid option MUST be linked to at least one attribute specified in the criteria
    # otherwise, there will be no data to display if that particular option is selected because
    # the intersection will be null.
    valid_options = []
    for v in candidate_options:
        loop_back_options = snapshot.filter_relationships.get(v, {})
        all_loop_back_options = []
        for dim in loop_back_options:
            all_loop_back_options.extend(loop_back_options[dim])

        if len(set(all_loop_back_options).intersection(all_criteria_attributes)) > 0:
            valid_options.append(v)

    # remove the prefix from each valid option and return the result
    return [option.split("__")[1] for option in valid_options]


def depluralize(x):
    return x[:-1] if x[-1] == "s" else x


def setup_retry_session(retries=3, backoff_factor=2, status_forcelist=(500, 502, 503, 504), method_whitelist=None):
    session = requests.Session()

    if method_whitelist is None:
        method_whitelist = {"GET"}

    retry = Retry(
        total=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
        allowed_methods=method_whitelist,
    )

    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)

    return session


def get_datasets_from_discover_api():
    # hardcode to staging backend if deployment is rdev or test
    deployment_stage = os.environ.get("DEPLOYMENT_STAGE")
    API_URL = DEPLOYMENT_STAGE_TO_API_URL.get(
        deployment_stage, "https://api.cellxgene.staging.single-cell.czi.technology"
    )

    datasets = {}
    if API_URL:
        session = setup_retry_session()
        dataset_metadata_url = f"{API_URL}/curation/v1/datasets?schema_version={CENSUS_CUBE_PINNED_SCHEMA_VERSION}"
        response = session.get(dataset_metadata_url)
        if response.status_code == 200:
            datasets = response.json()
    return datasets


def get_collections_from_discover_api():
    # hardcode to staging backend if deployment is rdev or test
    deployment_stage = os.environ.get("DEPLOYMENT_STAGE")
    API_URL = DEPLOYMENT_STAGE_TO_API_URL.get(
        deployment_stage, "https://api.cellxgene.staging.single-cell.czi.technology"
    )

    collections = {}
    if API_URL:
        session = setup_retry_session()
        dataset_metadata_url = f"{API_URL}/curation/v1/collections"
        response = session.get(dataset_metadata_url)
        if response.status_code == 200:
            collections = response.json()
    return collections


def build_filter_relationships(cell_counts_df: pd.DataFrame):
    # get a dataframe of the columns that are not numeric
    df_filters = cell_counts_df.select_dtypes(exclude="number")
    # get a numpy array of the column names with shape (1, n_cols)
    cols = df_filters.columns.values[None, :]

    # tile the column names row-wise to match the shape of the dataframe and concatenate to the values
    # this ensures that filter values will never collide across columns.
    mat = np.tile(cols, (cell_counts_df.shape[0], 1)) + "__" + df_filters.values

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


def to_dict(a, b):
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
    d = dict(zip(np.unique(a), [list(set(x)) for x in slists], strict=False))
    return d


def get_all_cell_type_ids_in_corpus(snapshot: CensusCubeSnapshot, root_node="CL:0000000") -> list[str]:
    """
    Retrieve all cell type IDs in the corpus that have at least one cell present, starting from a specified root node in the ontology.

    This function uses the ontology parser to get all descendant cell type IDs from the root node, including the root node itself.
    It then fetches the cell counts dataframe from the snapshot and groups it by cell type ontology term ID to sum up the cell counts.
    Cell types with zero cells are added to ensure all cell types from the ontology are represented in the output, even if they have no cells.
    Descendant cell counts are rolled up into each cell type node and remaining cell types with zero cells after rollup are filtered out.
    The function finally returns a list of cell type IDs that have cells.

    Parameters:
        snapshot (CensusCubeSnapshot): The snapshot to load.
        root_node (str): The ontology term ID of the root node from which to start gathering cell type IDs.

    Returns:
        list[str]: A list of cell type ontology term IDs that have at least one cell in the corpus.
    """

    all_cell_type_ids = ontology_parser.get_term_descendants(root_node, include_self=True)
    cell_counts_df = snapshot.cell_counts_df

    cell_counts_df = (
        cell_counts_df.groupby("cell_type_ontology_term_id").sum(numeric_only=True)[["n_cells"]].reset_index()
    )
    # to_attach is a DataFrame that will contain cell type ontology term ids that are not present in the input cell counts dataframe. These will be added with 0 counts.
    to_attach = pd.DataFrame()
    to_attach["cell_type_ontology_term_id"] = [
        i for i in all_cell_type_ids if i not in cell_counts_df["cell_type_ontology_term_id"].values
    ]
    to_attach["n_cells"] = 0
    cell_counts_df = pd.concat([cell_counts_df, to_attach], axis=0)
    cell_counts_df_rollup = rollup_across_cell_type_descendants(cell_counts_df)
    cell_counts_df_rollup = cell_counts_df_rollup.set_index("cell_type_ontology_term_id")["n_cells"]
    all_cell_type_ids_in_corpus = cell_counts_df_rollup.index.values[cell_counts_df_rollup.values > 0]
    return list(all_cell_type_ids_in_corpus)


def get_all_tissue_ids_in_corpus(snapshot: CensusCubeSnapshot) -> list[str]:
    """
    Retrieve all tissue IDs in the corpus that have at least one cell present.

    This function fetches the cell counts dataframe from the snapshot.
    It then uses a utility function to convert the grouped data into a dictionary mapping tissue IDs to cell type IDs.
    The function finally returns a list of tissue IDs that have cells.

    Parameters:
        snapshot (CensusCubeSnapshot): The snapshot to load.
    Returns:
        list[str]: A list of tissue ontology term IDs that have at least one cell in the corpus.
    """
    cell_counts_df = snapshot.cell_counts_df
    uberon_by_celltype = to_dict(
        cell_counts_df["tissue_ontology_term_id"], cell_counts_df["cell_type_ontology_term_id"]
    )
    return list(uberon_by_celltype.keys())


# cache finding descendants per cell type
@lru_cache(maxsize=None)
def descendants(cell_type):
    try:
        return ontology_parser.get_term_descendants(cell_type, include_self=True)
    except ValueError:
        return [cell_type]


@lru_cache(maxsize=None)
def ancestors(cell_type):
    try:
        return ontology_parser.get_term_ancestors(cell_type, include_self=True)
    except ValueError:
        return [cell_type]


def get_valid_descendants(
    cell_type: str, valid_cell_types: frozenset[str], cell_counts: Optional[dict[str, int]] = None
):
    """
    Get valid descendants for a cell type. Here, "validity" can mean one of two things:
    1. The descendant is present in the input list of valid cell types
    2. The descendant is present in the input list of valid cell types and has a different number of cells
    than the input cell type. Here, "number of cells" refers to the number of cells of that type AFTER roll-up.

    Terms can be suffixed with ";;{A_0}--{B_0}--..." to indicate that only descendants within the same
    group specified by the suffix should be considered. For example, if the term is
    "CL:0000540;;{tissue_0}--{disease_0}" then only descendants with the same tissue_0 and disease_0
    will be considered.

    If cell counts are provided, then a cell type has no valid descendants (even itself) if it has the same
    number of cells as one of its descendants. In these cases, the cell type is a redundant node and should
    be excluded. Redundant nodes are indicated by an empty list.

    Arguments
    ---------
    cell_type : str
        Cell type (cell type ontology term ID, potentially suffixed)
    valid_cell_types : frozenset[str]
        FrozenSet of valid cell types
    cell_counts : dict[str, int], optional, default=None
        Dictionary mapping cell type ontology term IDs (potentially suffixed) to the number of cells of that type.
        These cell counts are POST roll-up. This is used to exclude all descendants of a cell type if the cell type
        has the same number of cells as one of its descendants. In these cases, the cell type is a redundant node
        and should be excluded.


    Returns
    -------
    list[str] - An empty list indicates that the cell type is a redundant node and should be excluded.
    """
    # if the input cell types are suffixed, only consider descendants within the same group
    if ";;" in cell_type:
        prefix, suffix = cell_type.split(";;")
        suffix = ";;" + suffix
    else:
        prefix = cell_type
        suffix = ""

    # find descendants of cell type and re-suffix them
    relatives = [f"{i}{suffix}" for i in descendants(prefix)]
    valid_relatives = list(valid_cell_types.intersection(relatives))
    if cell_counts is not None:
        relative_cell_counts = [
            cell_counts.get(relative)
            for relative in valid_relatives
            if relative != cell_type and cell_counts.get(relative) is not None
        ]
        # if the cell type has the same number of cells as one of its descendants, then it is a redundant node
        # and has no valid relatives (including itself)
        if cell_counts.get(cell_type) in relative_cell_counts:
            valid_relatives = []

    return valid_relatives


def find_descendants_per_cell_type(cell_types):
    """
    Find the descendants for each cell type in the input list.
    Terms are only considered descendants if they are present in the input list.

    Terms can be suffixed with ";;{A_0}--{B_0}--..." to indicate that only descendants within the same
    group specified by the suffix should be considered. For example, if the term is
    "CL:0000540;;{tissue_0}--{disease_0}" then only descendants with the same tissue_0 and disease_0
    will be considered.

    This is useful for rolling up cell types within groups (e.g. tissues) but not across groups.


    Parameters
    ----------
    cell_types : list
        List of cell types (cell type ontology term IDs) to find descendants for

    Returns
    -------
    descendants_per_cell_type : list
        List of lists of descendants for each cell type in the input list.
    """
    lookup_table = {}
    cell_types_set = frozenset(cell_types)
    valid_descendants = []
    for cell_type in cell_types:
        if cell_type not in lookup_table:
            lookup_table[cell_type] = get_valid_descendants(cell_type, cell_types_set)
        valid_descendants.append(lookup_table[cell_type])
    return valid_descendants


def are_cell_types_not_redundant_nodes(cell_types, cell_counts):
    """
    Determines whether each cell type in a list of cell types is not a redundant node.
    Redundance means that a cell type has the same number of cells as one of its descendants. This indicates that
    the cell type contains exactly the same set of cells as one of its descendant subtrees. The only way this can occur is
    if the cell type belongs to a linear chain of cell types, all of which have the same number of cells. In these cases,
    the cell type is redundant and should be excluded. We only wish to keep the leaf nodes of these linear chains.

    Args:
    - cell_types (list of str): A list of cell type names.
    - cell_counts (dict): A dictionary mapping cell type names to the number of cells of that type.

    Returns:
    - is_not_redundant (list of bool): A list of boolean values indicating whether each cell type is a redundant node.
    """
    lookup_table = {}
    is_not_redundant = []
    cell_types_set = frozenset(cell_types)
    for cell_type in cell_types:
        if cell_type not in lookup_table:
            relatives = get_valid_descendants(cell_type, cell_types_set, cell_counts=cell_counts)
            lookup_table[cell_type] = len(relatives) > 0
        is_not_redundant.append(lookup_table[cell_type])
    return is_not_redundant


def are_cell_types_colinear(cell_type1, cell_type2):
    """
    Determine if two cell types are colinear in the ontology.
    Colinearity means that cell type 1 is an ancestor of cell type 2
    or vice-versa.
    Arguments
    ---------
    cell_type1 : str
        Cell type 1 (cell type ontology term id)
    cell_type2 : str
        Cell type 2 (cell type ontology term id)
    Returns
    -------
    bool
    """
    descendants1 = descendants(cell_type1)
    descendants2 = descendants(cell_type2)
    ancestors1 = ancestors(cell_type1)
    ancestors2 = ancestors(cell_type2)
    return len(set(descendants1).intersection(ancestors2)) > 0 or len(set(descendants2).intersection(ancestors1)) > 0


def get_overlapping_cell_type_descendants(cell_type1, cell_type2):
    """
    Get overlapping cell type descendants

    Arguments
    ---------
    cell_type1 : str
        Cell type 1 (cell type ontology term id)
    cell_type2 : str
        Cell type 2 (cell type ontology term id)
    Returns
    -------
    list[str]
    """
    descendants1 = descendants(cell_type1)
    descendants2 = descendants(cell_type2)

    return list(set(descendants1).intersection(descendants2))


def rollup_across_cell_type_descendants(
    df, cell_type_col="cell_type_ontology_term_id", parallel=True, ignore_cols=None
) -> pd.DataFrame:
    """
    Aggregate values for each cell type across its descendants in the input dataframe.

    The non-numeric columns in the input dataframe must contain cell type ontology term IDs,
    and are treated as the dimensions of a multi-dimensional numpy array. The numeric data in
    the dataframe is slotted into this array and rolled up along the first axis (which will always
    correspond to the cell type ontology term IDs). The resulting rolled up array is reshaped back
    into the tidy dataframe and returned.

    This ensures that cell types are only rolled up within each combination of other dimensions
    (e.g. tissue, gene, organism). We wouldn't want to roll up expressions across genes and tissues.

    Parameters
    ----------
    df : pandas DataFrame
        Tidy dataframe containing the dimensions across which the numeric columns will be
        aggregated. The dataframe must have a column containing the cell type ontology term IDs.
        By default, the column name is "cell_type_ontology_term_id".

    cell_type_col : str, optional, default="cell_type_ontology_term_id"
        Name of the column in the input dataframe containing the cell type ontology term IDs.

    parallel : bool, optional, default=True
        If True, uses numba's `prange` to parallelize the rollup operation.
        Set to False if invoking this function in parallel subprocesses.

    ignore_cols : list, optional, default=None
        List of column names to ignore when rolling up the numeric columns.

    Returns
    -------
    df : pandas DataFrame
        Tidy dataframe with the same dimensions as the input dataframe, but with the numeric
        columns aggregated across the cell type's descendants.
    """
    df = df.copy()
    # numeric data
    numeric_df = df.select_dtypes(include="number")
    # non-numeric data
    dimensions_df = df.select_dtypes(exclude="number")
    # move the cell type column to the front of the dataframe
    cell_type_column = dimensions_df.pop(cell_type_col)
    dimensions_df.insert(0, cell_type_col, cell_type_column)

    # calculate integer indices for each non-numeric column in the input dataframe
    # and calculate the shape of the output array
    dim_indices = []
    dim_shapes = []
    for col in dimensions_df.columns:
        # map each unique value in the column to an integer index
        dim = dimensions_df[col].to_numpy()
        unique_dim = dimensions_df[col].unique()
        indices = pd.Series(data=range(len(unique_dim)), index=unique_dim)[dim].to_numpy()
        dim_shapes.append(len(unique_dim))
        dim_indices.append(indices)
    # the last dimension corresponds to the numeric columns
    dim_shapes.append(numeric_df.shape[1])
    dim_shapes = tuple(dim_shapes)
    array_to_sum = np.zeros(dim_shapes)
    # slot the numeric data into the multi-dimensional numpy array

    array_to_sum[tuple(dim_indices)] = numeric_df.to_numpy()
    cell_types = cell_type_column.unique()

    summed = rollup_across_cell_type_descendants_array(array_to_sum, cell_types, parallel=parallel)
    # extract numeric data and write back into the dataframe
    summed = summed[tuple(dim_indices)]
    dtypes = numeric_df.dtypes
    for col, array in zip(numeric_df.columns, summed.T, strict=False):
        if ignore_cols and col not in ignore_cols or not ignore_cols:
            df[col] = array.astype(dtypes[col])

    return df


def rollup_across_cell_type_descendants_array(array_to_sum, cell_types, parallel=True) -> np.ndarray:
    """
    Aggregate values for each cell type across its descendants in the input array.
    Cell types must be the first dimension of the input array.

    Parameters
    ----------
    array_to_sum : numpy array
        Multi-dimensional numpy array containing the numeric data to be rolled up. The first
        dimension must correspond to the cell type ontology term IDs.

    cell_types : list
        List of cell type ontology term IDs corresponding to the first dimension of the input

    parallel : bool, optional, default=True
        If True, uses numba's `prange` to parallelize the rollup operation.
        Set to False if invoking this function in parallel subprocesses.

    Returns
    -------
    summed : numpy array
        Multi-dimensional numpy array aggregated across the cell type's descendants.
    """
    descendants = find_descendants_per_cell_type(cell_types)
    # a pandas series to map cell types to their index in the input arrays
    indexer = pd.Series(index=cell_types, data=range(len(cell_types)))
    descendants_indexes = [indexer[children].to_numpy() for children in descendants]
    # flatten the descendant indices into a single array and create a linear
    # index array for slicing out the descendants per cell type. The array
    # must be flattened to satisfy numba type requirements.
    z = 0
    linear_indices = [z]
    for ix in descendants_indexes:
        z += len(ix)
        linear_indices.append(z)
    linear_indices = np.array(linear_indices)
    descendants_indexes = np.concatenate(descendants_indexes)

    # roll up the multi-dimensional array across cell types (first axis)
    summed = np.zeros_like(array_to_sum)
    if parallel:
        _sum_array_elements__parallel(array_to_sum, summed, descendants_indexes, linear_indices)
    else:
        _sum_array_elements(array_to_sum, summed, descendants_indexes, linear_indices)
    return summed


@nb.njit(parallel=True, fastmath=True, nogil=True)
def _sum_array_elements__parallel(array, summed, descendants_indexes, linear_indices):
    for i in nb.prange(len(linear_indices) - 1):
        index = descendants_indexes[linear_indices[i] : linear_indices[i + 1]]
        for j in index:
            summed[i] += array[j]


@nb.njit(fastmath=True, nogil=True)
def _sum_array_elements(array, summed, descendants_indexes, linear_indices):
    for i in range(len(linear_indices) - 1):
        index = descendants_indexes[linear_indices[i] : linear_indices[i + 1]]
        for j in index:
            summed[i] += array[j]
