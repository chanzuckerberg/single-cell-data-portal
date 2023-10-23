from functools import lru_cache

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


def find_descendants_per_cell_type(cell_types):
    """
    Find the ancestors or descendants for each cell type in the input list.

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
    cell_types_set = set(cell_types)
    relatives_per_cell_type = []
    for cell_type in cell_types:
        if ";;" in cell_type:
            prefix, suffix = cell_type.split(";;")
            suffix = ";;" + suffix
        else:
            prefix = cell_type
            suffix = ""
        if cell_type not in lookup_table:
            relatives = [f"{i}{suffix}" for i in descendants(prefix)]
            lookup_table[cell_type] = list(cell_types_set.intersection(relatives))

        relatives_per_cell_type.append(lookup_table[cell_type])

    return relatives_per_cell_type


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
