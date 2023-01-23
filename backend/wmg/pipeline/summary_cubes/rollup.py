import owlready2
import pandas as pd
import numpy as np
from backend.wmg.data.constants import CL_BASIC_PERMANENT_URL_OWL


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
    onto = owlready2.get_ontology(CL_BASIC_PERMANENT_URL_OWL)
    onto.load()

    descendants_per_cell_type = []
    for cell_type in cell_types:
        cell_type_iri = cell_type.replace(":", "_")
        entity = onto.search_one(iri=f"http://purl.obolibrary.org/obo/{cell_type_iri}")
        if entity:
            descendants = [i.name.replace("_", ":") for i in entity.descendants()]
        else:
            descendants = [cell_type]
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
    summed_arrays = []
    for array in arrays_to_sum:
        summed = np.zeros_like(array)
        for i, children in enumerate(descendants):
            # indexer is used to map children cell types to their index
            # in the input arrays
            summed[i] += array[indexer[children]].sum(axis=0)
        summed_arrays.append(summed)
    return summed_arrays
