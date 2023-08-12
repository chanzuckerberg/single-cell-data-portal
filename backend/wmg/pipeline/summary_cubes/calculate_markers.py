import json
from functools import lru_cache, wraps
from typing import Optional, Union

import numpy as np
import pandas as pd
import tiledb
from numba import njit, prange
from scipy import stats

from backend.common.utils.exceptions import MarkerGeneCalculationException
from backend.wmg.data.query import FmgQueryCriteria, WmgCubeQueryParams, WmgQuery
from backend.wmg.data.rollup import (
    are_cell_types_colinear,
    descendants,
    rollup_across_cell_type_descendants,
    rollup_across_cell_type_descendants_array,
)
from backend.wmg.data.schemas.data_schema_config import (
    WMG_DATA_SCHEMA_VERSION,
    WRITER_WMG_CUBE_QUERY_VALID_ATTRIBUTES,
    WRITER_WMG_CUBE_QUERY_VALID_DIMENSIONS,
)
from backend.wmg.data.snapshot import (
    CELL_COUNTS_CUBE_NAME,
    DATASET_TO_GENE_IDS_FILENAME,
    EXPRESSION_SUMMARY_FMG_CUBE_NAME,
    WmgSnapshot,
    load_snapshot,
)


@njit(parallel=True)
def nanpercentile_2d(arr, percentile, axis):
    """
    Calculate the specified percentile of a 2D array along an axis, ignoring NaN values.

    Parameters:
        arr: 2D array to calculate percentile of
        percentile: percentile to calculate, as a number between 0 and 100
        axis: axis along which to calculate percentile

    Returns:
        The specified percentile of the 2D array along the specified axis.
    """
    if axis == 0:
        result = np.empty(arr.shape[1])
        for i in prange(arr.shape[1]):
            arr_column = arr[:, i]
            result[i] = nanpercentile(arr_column, percentile)
        return result
    else:
        result = np.empty(arr.shape[0])
        for i in prange(arr.shape[0]):
            arr_row = arr[i, :]
            result[i] = nanpercentile(arr_row, percentile)
        return result


@njit
def nanpercentile(arr, percentile):
    arr_without_nan = arr[np.logical_not(np.isnan(arr))]
    length = len(arr_without_nan)

    if length == 0:
        return np.nan

    return np.percentile(arr_without_nan, percentile)


def _make_hashable(func):
    """
    Implicitly convert unhashable data structures (list and dict) to hashable strings
    for memoization.

    WmgSnapshot objects are not hashable, so we ignore them for memoization.
    """

    class HDict(dict):
        def __hash__(self):
            return hash(json.dumps(dict(sorted(self.items()))))

    class HList(list):
        def __hash__(self):
            return hash(json.dumps(self))

    @wraps(func)
    def wrapped(*args, **kwargs):
        new_args = []
        for arg in args:
            if isinstance(arg, dict):
                arg = HDict(arg)
            elif isinstance(arg, list):
                arg = HList(arg)
            new_args.append(arg)

        new_kwargs = {}
        for k, v in kwargs.items():
            if isinstance(v, dict):
                v = HDict(v)
            elif isinstance(v, list):
                v = HList(v)
            new_kwargs[k] = v

        return func(*new_args, **new_kwargs)

    return wrapped

    # this decorator function is required to allow lru_cache to


# memoize the function as dicts, lists, and WmgSnapshot objects
# are not typically hashable.
@_make_hashable
@lru_cache(maxsize=None)
def _query_tiledb_context(
    cell_type: str,
    tissue: str,
    organism: str,
    corpus: Optional[Union[WmgSnapshot, str]] = None,
):
    if isinstance(corpus, str):
        expression_summary_fmg_cube = tiledb.open(f"{corpus}/{EXPRESSION_SUMMARY_FMG_CUBE_NAME}")
        cell_counts_cube = tiledb.open(f"{corpus}/{CELL_COUNTS_CUBE_NAME}")
        with open(f"{corpus}/{DATASET_TO_GENE_IDS_FILENAME}") as fp:
            dataset_to_gene_ids = json.load(fp)
    else:
        if corpus is None:
            corpus = load_snapshot(
                snapshot_schema_version=WMG_DATA_SCHEMA_VERSION,
                read_versioned_snapshot=True,
            )

        assert isinstance(corpus, WmgSnapshot)
        expression_summary_fmg_cube = corpus.expression_summary_fmg_cube
        cell_counts_cube = corpus.cell_counts_cube
        dataset_to_gene_ids = corpus.dataset_to_gene_ids

    # build a snapshot containing only the relevant cubes with everything else set to None.
    # this is to use the WmgQuery interface.
    snapshot = WmgSnapshot(
        snapshot_identifier="",
        expression_summary_fmg_cube=expression_summary_fmg_cube,
        cell_counts_cube=cell_counts_cube,
        dataset_to_gene_ids=dataset_to_gene_ids,
        marker_genes_cube=None,
        expression_summary_cube=None,
        expression_summary_default_cube=None,
        cell_type_orderings=None,
        primary_filter_dimensions=None,
        filter_relationships=None,
    )
    criteria = FmgQueryCriteria(**{"tissue_ontology_term_ids": [tissue], "organism_ontology_term_id": organism})

    cube_query_params = WmgCubeQueryParams(
        cube_query_valid_attrs=WRITER_WMG_CUBE_QUERY_VALID_ATTRIBUTES,
        cube_query_valid_dims=WRITER_WMG_CUBE_QUERY_VALID_DIMENSIONS,
    )

    q = WmgQuery(snapshot, cube_query_params)
    query = q.expression_summary_fmg(criteria)

    cell_counts_query = q.cell_counts(criteria)

    criteria.cell_type_ontology_term_ids = descendants(cell_type)
    query = q.expression_summary_fmg(criteria)
    # group-by and sum
    agg = query.groupby(["cell_type_ontology_term_id", "gene_ontology_term_id"]).sum(numeric_only=True).reset_index()
    agg = rollup_across_cell_type_descendants(agg)
    # remove tidy
    agg = agg.groupby(["cell_type_ontology_term_id", "gene_ontology_term_id"]).sum(numeric_only=True)

    n_cells = cell_counts_query.groupby(["dataset_id", "cell_type_ontology_term_id"]).sum(numeric_only=True)[
        "n_total_cells"
    ]
    genes = list(agg.index.levels[1])
    n_cells_per_gene, n_cells_index = _calculate_true_n_cells(n_cells, genes, dataset_to_gene_ids)

    return agg, n_cells_per_gene, n_cells_index, genes


def _query_target(
    cell_type: str,
    context_agg: pd.DataFrame,
    n_cells_per_gene_context: np.array,
    n_cells_index_context: pd.Index,
):
    """
    Extract data for the target population from the context query results.

    Arguments
    ---------
    cell_type - str
        Cell type of the target population

    context_agg - pd.DataFrame
        Dataframe of aggregated metrics from the context query

    n_cells_per_gene_context - numpy array
        Each row in the array contains the true number of cells present for the combination of
        filters in the context query.

    n_cells_index_context - pandas Index
        Each element in the index contains the cell type.

    Returns
    -------
    - agg: target query result aggregated by filter keys
    - n_cells_per_gene_target: numpy array
        each row in the array contains the true number of cells present for the target population.
    """

    # find mismatch between target_filter and context_filter
    # comparisons will be made across the mismatched dimensions

    filt = context_agg.index.get_level_values("cell_type_ontology_term_id") == cell_type
    filt_ncells = n_cells_index_context == cell_type

    if not np.any(filt) or not np.any(filt_ncells):
        raise MarkerGeneCalculationException("No cells found for target population.")

    target_agg = context_agg[filt]
    n_cells_per_gene_target = n_cells_per_gene_context[filt_ncells].sum(axis=0, keepdims=True)
    return target_agg, n_cells_per_gene_target


def _calculate_true_n_cells(n_cells, genes, dataset_to_gene_ids):
    """
    Calculates the true number of cells per gene for each combination of filters.

    Arguments
    ---------
    n_cells - pandas.Series
        A pandas series with the number of cells per combination of filters.
        The index of the series is a multi-index with the filter keys as levels.

    genes - list
        A list of genes to calculate true population sizes for.

    dataset_to_gene_ids - dict
        A dictionary mapping dataset IDs to a list of gene IDs.

    Returns
    -------
    n_cells_array - numpy.ndarray
        A numpy array with the number of cells per gene for each combination of filters.

    index - pandas.MultiIndex
        A pandas multi-index with the filter keys as levels corresponding to each row of
        the returned array.
    """
    genes_indexer = pd.Series(index=genes, data=range(len(genes)))

    # fill in an array stratified across desired filter combinations + dataset ID
    # with the number of cells in each group for genes that are PRESENT in the dataset
    t_n_cells = np.zeros((n_cells.shape[0], len(genes)))
    for i, index in enumerate(n_cells.index):
        n = n_cells[index]
        dataset_id = index[0]
        present = list(set(dataset_to_gene_ids[dataset_id]).intersection(genes))
        t_n_cells[i, genes_indexer[present]] = n

    # get the tuples of filter values, excluding dataset IDs if they are not specified
    # as a filter by the user.
    cell_types = n_cells.index.get_level_values("cell_type_ontology_term_id")

    # sum up the cell count arrays across duplicate groups
    # (groups can be duplicate after excluding dataset_ids)
    t_n_cells_sum = {}
    for i, k in enumerate(cell_types):
        summer = t_n_cells_sum.get(k, np.zeros(t_n_cells.shape[1]))
        summer += t_n_cells[i]
        t_n_cells_sum[k] = summer

    n_cells_array = np.vstack(list(t_n_cells_sum.values()))

    n_cells_array = rollup_across_cell_type_descendants_array(n_cells_array, cell_types)
    index = pd.Index(cell_types, name="cell_type_ontology_term_id")
    return n_cells_array, index


def _run_ttest(sum1, sumsq1, n1, sum2, sumsq2, n2):
    """
    Calculate t-test statistics.

    Arguments
    ---------
    sum1 - np.ndarray (1 x M)
        Array of sum expression in target pop for each gene
    sumsq1 - np.ndarray (1 x M)
        Array of sum of squared expressions in target pop for each gene
    n1 - np.ndarray (1 x M)
        Number of cells in target pop for each gene
    sum2 - np.ndarray (N x M)
        Array of sum expression in context pops for each gene
    sumsq2 - np.ndarray (N x M)
        Array of sum of squared expressions in context pops for each gene
    n2 - np.ndarray (N x M)
        Number of cells in context pops for each gene

    Returns
    -------
    pvals_adj - np.ndarray
        adjusted p-values for each comparison for each gene
    effects - np.ndarray
        effect sizes for each comparison for each gene
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        mean1 = sum1 / n1
        meansq1 = sumsq1 / n1

        mean2 = sum2 / n2
        meansq2 = sumsq2 / n2

        var1 = meansq1 - mean1**2
        var1[var1 < 0] = 0
        var2 = meansq2 - mean2**2
        var2[var2 < 0] = 0

        var1_n = var1 / n1
        var2_n = var2 / n2
        sum_var_n = var1_n + var2_n
        dof = sum_var_n**2 / (var1_n**2 / (n1 - 1) + var2_n**2 / (n2 - 1))
        tscores = (mean1 - mean2) / np.sqrt(sum_var_n)
        effects = (mean1 - mean2) / np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 1))

    pvals = stats.t.sf(tscores, dof)
    return pvals, effects


def _post_process_stats(
    cell_type_target,
    cell_types_context,
    genes,
    pvals,
    effects,
    nnz,
    test="ttest",
    min_num_expr_cells=25,
    percentile=0.15,
    n_markers=10,
):
    """
    Process and aggregate statistics into the output dictionary format.

    Arguments
    ---------
    cell_type_target - str or None
        Cell type ontology term ID corresponding to the target population. If
        None, a cell type was not specified for the target population and so marker
        genes are not being calculated for cell types.

    cell_types_context - list or None
        List of cell types in the context populations. If None, marker genes are not
        being calculated for cell types.

    genes - list
        List of genes corresponding to each p-value and effect size

    pvals - np.ndarray
        N x M array of p-values for each comparison and each gene

    effects - np.ndarray
        N x M array of effect sizes for each comparison and each gene

    nnz - np.ndarray
        1 x M array of number of nonzero expressions for target population (> 0)

    test - str, optional, default "ttest"
        The statistical test used. Historically "ttest" or "binomtest" were both supported, currently
        only "ttest" is supported.

    min_num_expr_cells - int, optional, default 25
        The minimum number of nonzero expressing cells required for marker genes

    percentile - float, optional, default 0.15
        The percentile of effect sizes to select as the representative effect size.

    n_markers - int, optional, default 10
        Number of top markers to return. If None, all marker genes with effect size > 0 are returned.
    """
    zero_out = nnz.flatten() < min_num_expr_cells
    effects[:, zero_out] = 0
    pvals[:, zero_out] = 1

    if cell_type_target is not None and cell_types_context is not None:
        is_colinear = np.array(
            [are_cell_types_colinear(cell_type, cell_type_target) for cell_type in cell_types_context]
        )
        effects[is_colinear] = np.nan
        pvals[is_colinear] = np.nan
        pvals[:, np.all(np.isnan(effects), axis=0)] = 1
        effects[:, np.all(np.isnan(effects), axis=0)] = 0

    # aggregate
    effects = nanpercentile_2d(effects, percentile * 100, 0)
    pvals = np.sort(pvals, axis=0)[int(np.round(0.05 * pvals.shape[0]))]
    effects[zero_out] = np.nan
    pvals[zero_out] = np.nan

    if n_markers:
        markers = np.array(genes)[np.argsort(-effects)[:n_markers]]
        p = pvals[np.argsort(-effects)[:n_markers]]
        effects = effects[np.argsort(-effects)[:n_markers]]
    else:
        markers = np.array(genes)[np.argsort(-effects)]
        p = pvals[np.argsort(-effects)]
        effects = effects[np.argsort(-effects)]
    statistics = []
    final_markers = []
    for i in range(len(p)):
        pi = p[i]
        ei = effects[i]
        if ei is not np.nan and pi is not np.nan:
            statistics.append({f"p_value_{test}": pi, f"effect_size_{test}": ei})
            final_markers.append(markers[i])
    return dict(zip(list(final_markers), statistics))


def get_markers(cell_type, tissue, organism, corpus=None, n_markers=10, percentile=0.15):
    """
    Calculate marker genes using the t-test.

    Arguments
    ---------
    cell_type - str
        cell type of interest

    tissue - dict
        tissue of interest

    corpus - str or WmgSnapshot, optional, default None
        If string, it is the path to the snapshot.
        If WmgSnapshot, it is the snapshot object.
        If None, the snapshot will be fetched from AWS.

    n_markers - int, optional, default 10
        Number of top markers to return. If None, all marker genes with effect size > 0 are returned.

    percentile - float, optional, default 0.15
        The percentile of effect sizes to select as the representative effect size.

    Returns
    -------
    Dictionary of marker genes as keys and a dictionary of p-value and effect size as values.
    """

    """
    Calculate marker genes using the t-test.

    Arguments
    ---------
    cell_type - str
        cell type ontology term ID corresponding to the target population

    tissue - str
        tissue ontology term ID corresponding to the target population

    corpus - str or WmgSnapshot, optional, default None
        If string, it is the path to the snapshot.
        If WmgSnapshot, it is the snapshot object.
        If None, the snapshot will be fetched from AWS.

    n_markers - int, optional, default 10
        Number of top markers to return. If None, all marker genes with effect size > 0 are returned.

    percentile - float, optional, default 0.15
        The percentile of effect sizes to select as the representative effect size.

    Returns
    -------
    Dictionary of marker genes as keys and a dictionary of p-value and effect size as values.
    """

    context_agg, n_context, cell_types_index_context, genes = _query_tiledb_context(
        cell_type, tissue, organism, corpus=corpus
    )
    target_agg, n_target = _query_target(cell_type, context_agg, n_context, cell_types_index_context)

    target_agg = target_agg.groupby("gene_ontology_term_id").sum(numeric_only=True)

    genes_target = list(target_agg.index)

    genes_context = list(context_agg.index.get_level_values("gene_ontology_term_id"))
    cell_types_context = list(context_agg.index.get_level_values("cell_type_ontology_term_id"))
    cell_types_indexer = pd.Series(index=cell_types_index_context, data=range(len(cell_types_index_context)))
    genes_indexer = pd.Series(index=genes, data=np.arange(len(genes)))

    cell_types_indices_context = list(cell_types_indexer[cell_types_context])
    genes_indices_context = list(genes_indexer[genes_context])
    genes_indices_target = list(genes_indexer[genes_target])

    target_data_nnz = np.zeros((1, len(genes)))
    target_data_nnz[0, genes_indices_target] = list(target_agg["nnz"])

    target_data_sum = np.zeros((1, len(genes)))
    target_data_sum[0, genes_indices_target] = list(target_agg["sum"])
    target_data_sumsq = np.zeros((1, len(genes)))
    target_data_sumsq[0, genes_indices_target] = list(target_agg["sqsum"])

    context_data_sum = np.zeros((len(cell_types_index_context), len(genes)))
    context_data_sum[cell_types_indices_context, genes_indices_context] = list(context_agg["sum"])
    context_data_sumsq = np.zeros((len(cell_types_index_context), len(genes)))
    context_data_sumsq[cell_types_indices_context, genes_indices_context] = list(context_agg["sqsum"])

    pvals, effects = _run_ttest(
        target_data_sum, target_data_sumsq, n_target, context_data_sum, context_data_sumsq, n_context
    )

    cell_types_context = cell_types_index_context.get_level_values("cell_type_ontology_term_id")

    return _post_process_stats(
        cell_type,
        cell_types_context,
        genes,
        pvals,
        effects,
        target_data_nnz,
        test="ttest",
        n_markers=n_markers,
        percentile=percentile,
    )
