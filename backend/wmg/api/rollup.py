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
    return descendants_per_cell_type


def rollup_across_cell_type_descendants(cell_types, arrays_to_sum):
    """
    Aggregate values for each cell type across its descendants in the input arrays.

    Parameters
    ----------
    cell_types : list
        List of cell types corresponding to each row of the input arrays

    arrays_to_sum : list
        List of numpy arrays to roll up across descendants. Each array
        should have the same length as cell_types.

    Returns
    -------
    summed_arrays : list
        List of numpy arrays with the same shape as the input arrays, but
        with values rolled up across descendants for each cell type (row).
    """

    for array in arrays_to_sum:
        assert len(array) == len(cell_types), "Input arrays must have the same length as cell_types"

    descendants = find_descendants_per_cell_type(cell_types)
    for i, children in enumerate(descendants):
        # remove descendent cell types that are not cell types in the WMG data
        descendants[i] = list(set(children).intersection(cell_types))

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

    summed_arrays = []
    for array in arrays_to_sum:
        summed = np.zeros_like(array)
        _array_summer(array, summed, descendants_indexes, linear_indices)
        summed_arrays.append(summed)
    return summed_arrays


@nb.njit(parallel=True, fastmath=True, nogil=True)
def _array_summer(array, summed, descendants_indexes, linear_indices):
    for i in nb.prange(len(linear_indices) - 1):
        index = descendants_indexes[linear_indices[i] : linear_indices[i + 1]]
        for j in index:
            summed[i] += array[j]
