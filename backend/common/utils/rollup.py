from functools import lru_cache
from typing import Optional

import numba as nb
import numpy as np
import owlready2
import pandas as pd

from backend.wmg.data.constants import CL_BASIC_OWL_NAME
from backend.wmg.data.utils import get_pinned_ontology_url

# ontology object
ontology = owlready2.get_ontology(get_pinned_ontology_url(CL_BASIC_OWL_NAME))
ontology.load()


# cache finding descendants per cell type
@lru_cache(maxsize=None)
def descendants(cell_type):
    global ontology
    cell_type_iri = cell_type.replace(":", "_")
    entity = ontology.search_one(iri=f"http://purl.obolibrary.org/obo/{cell_type_iri}")
    descendants = [i.name.replace("_", ":") for i in entity.descendants()] if entity else [cell_type]
    return descendants


@lru_cache(maxsize=None)
def ancestors(cell_type):
    global ontology
    cell_type_iri = cell_type.replace(":", "_")
    entity = ontology.search_one(iri=f"http://purl.obolibrary.org/obo/{cell_type_iri}")
    ancestors = [i.name.replace("_", ":") for i in entity.ancestors()] if entity else [cell_type]
    return ancestors


def get_valid_descendants(
    cell_type: str, valid_cell_types: frozenset[str], cell_counts: Optional[dict[str, int]] = None
):
    """
    Get valid descendants for a cell type. Here, "validity" can mean one of two things:
    1. The descendant is present in the input list of valid cell types
    2. The descendant is present in the input list of valid cell types and has a different number of cells
    than the input cell type.

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
        This is used to exclude all descendants of a cell type if the cell type has the same number of cells as one
        of its descendants. In these cases, the cell type is a redundant node and should be excluded.


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
    Colinearity means that cell type 1 is an aacestor of cell type 2
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
    for col, array in zip(numeric_df.columns, summed.T):
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
