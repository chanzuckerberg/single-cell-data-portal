from pronto import Ontology
import pandas as pd
import numpy as np
from backend.wmg.data.constants import CL_BASIC_PERMANENT_URL


def aggregate_across_cell_type_descendants(cell_types, arrays_to_sum):
    """
    Aggregate values for each cell type across its descendants in the input arrays.

    Parameters
    ----------
    cell_types : list
        List of cell types corresponding to each row of the input arrays

    arrays_to_sum : list
        List of numpy arrays to aggregate across descendants. Each array
        should have the same number of rows as cell_types.

    Returns
    -------
    summed_arrays : list
        List of numpy arrays with the same shape as the input arrays, but
        with values aggregated across descendants for each cell type (row).
    """
    onto = Ontology(CL_BASIC_PERMANENT_URL)
    descendants_for_node = lambda cl: list(onto[cl].subclasses().to_set().ids)
    descendants = [descendants_for_node(i) for i in cell_types]
    for i, children in enumerate(descendants):
        descendants[i] = list(set(children).intersection(cell_types))
    indexer = pd.Series(index=cell_types, data=range(len(cell_types)))

    summed_arrays = []
    for array in arrays_to_sum:
        summed = np.zeros_like(array)
        flag = summed.ndim == 1
        if flag:
            summed = summed[:, None]
        for i, children in enumerate(descendants):
            summed[i] += array[indexer[children]].sum(0)
        if flag:
            summed = summed.flatten()
        summed_arrays.append(summed)
    return summed_arrays
