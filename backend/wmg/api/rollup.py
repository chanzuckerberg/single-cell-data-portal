import owlready2
import pandas as pd
import numpy as np
import numba as nb
from functools import lru_cache
from backend.wmg.data.constants import CL_BASIC_PERMANENT_URL_OWL


# ontology object
onto = None


def find_descendants_per_cell_type(cell_types):
    """
    Find descendants for each cell type in the input list.

    Parameters
    ----------
    cell_types : list
        List of cell types (cell type ontology term IDs) to find descendants for

    Returns
    -------
    descendants_per_cell_type : list
        List of lists of descendants for each cell type in the input list.
    """

    global onto
    if not onto:
        onto = owlready2.get_ontology(CL_BASIC_PERMANENT_URL_OWL)
        onto.load()

    # cache finding descendants per cell type
    @lru_cache(maxsize=None)
    def _descendants(cell_type):
        cell_type_iri = cell_type.replace(":", "_")
        entity = onto.search_one(iri=f"http://purl.obolibrary.org/obo/{cell_type_iri}")
        if entity:
            descendants = [i.name.replace("_", ":") for i in entity.descendants()]
        else:
            descendants = [cell_type]
        return descendants

    descendants_per_cell_type = []
    for cell_type in cell_types:
        descendants = _descendants(cell_type)
        descendants_per_cell_type.append(descendants)

    for i, children in enumerate(descendants_per_cell_type):
        # remove descendent cell types that are not cell types in the WMG data
        descendants_per_cell_type[i] = list(set(children).intersection(cell_types))

    return descendants_per_cell_type


def rollup_across_cell_type_descendants(df, cell_type_col="cell_type_ontology_term_id"):
    """
    Aggregate values for each cell type across its descendants in the input arrays.

    Parameters
    ----------
    df : pandas DataFrame
        Tidy dataframe containing the dimensions across which the numeric columns will be
        aggregated. The dataframe must have a column containing the cell type ontology term IDs.
        By default, the column name is "cell_type_ontology_term_id".

    cell_type_col : str, optional, default="cell_type_ontology_term_id"
        Name of the column in the input dataframe containing the cell type ontology term IDs.

    Returns
    -------
    df : pandas DataFrame
        Tidy dataframe with the same dimensions as the input dataframe, but with the numeric
        columns aggregated across the cell type's descendants.
    """
    df = df.copy()

    numeric_df = df.select_dtypes(include="number")
    dimensions_df = df.select_dtypes(exclude="number")

    cell_type_column = dimensions_df.pop(cell_type_col)
    dimensions_df.insert(0, cell_type_col, cell_type_column)

    dim_indices = []
    dim_shapes = []
    for col in dimensions_df.columns:
        dim = dimensions_df[col].to_numpy()
        unique_dim = dimensions_df[col].unique()
        indices = pd.Series(data=range(len(unique_dim)), index=unique_dim)[dim].to_numpy()
        dim_shapes.append(len(unique_dim))
        dim_indices.append(indices)
    dim_shapes.append(numeric_df.shape[1])
    dim_shapes = tuple(dim_shapes)
    array_to_sum = np.zeros(dim_shapes)
    array_to_sum[tuple(dim_indices)] = numeric_df.to_numpy()

    cell_types = cell_type_column.unique()

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

    summed = np.zeros_like(array_to_sum)
    _array_summer(array_to_sum, summed, descendants_indexes, linear_indices)

    summed = summed[tuple(dim_indices)]
    dtypes = numeric_df.dtypes
    for col, array in zip(numeric_df.columns, summed.T):
        df[col] = array.astype(dtypes[col])

    return df


@nb.njit(parallel=True, fastmath=True, nogil=True)
def _array_summer(array, summed, descendants_indexes, linear_indices):
    for i in nb.prange(len(linear_indices) - 1):
        index = descendants_indexes[linear_indices[i] : linear_indices[i + 1]]
        for j in index:
            summed[i] += array[j]
