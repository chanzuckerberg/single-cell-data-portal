import pandas as pd
import numpy as np
from scipy import stats
import json
import tiledb
from functools import lru_cache, wraps
from backend.wmg.data.query import FmgQueryCriteria, WmgQuery
from backend.wmg.data.snapshot import (
    load_snapshot,
    WmgSnapshot,
    EXPRESSION_SUMMARY_FMG_CUBE_NAME,
    CELL_COUNTS_CUBE_NAME,
    DATASET_TO_GENE_IDS_FILENAME,
)
from typing import Union, Optional


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
def _query_tiledb(
    filters: dict,
    corpus: Optional[Union[WmgSnapshot, str]] = None,
    group_by_dims: Optional[list] = None,
    genes: Optional[list] = None,
):
    """
    Query the tileDB array and return required data artifacts for downstream processing.

    Arguments
    ---------
    filters - dict,
        Dictionary of filters for the target population

    corpus - str or WmgSnapshot, optional, default None
        If string, it is the path to the snapshot.
        If WmgSnapshot, it is the snapshot object.
        If None, the snapshot will be fetched from AWS.

    group_by_dims - list, optional, default None
        A list of filter dimensions to aggregate the tileDB query result by.
        If group-by dimensions not provided, use the keys specified in the `filters`

    genes - list, optional, default None
        A list of genes to calculate true population sizes for.
        If None, the genes in the aggregated query results are used.

    Returns
    -------
    - agg: query result aggregated by filter keys
    - t_n_cells_sum: a dictionary of 1D numpy arrays with length # of genes
        each value is the true number of cells present for the combination of
        filters (keys).
        Ex: {(cell_type1,tissue1,organism1): np.array([0,100,100,0,80])}
    - genes: a list of all genes with nonzero expression found in the query result
    """

    criteria = FmgQueryCriteria(**filters)

    if isinstance(corpus, str):
        expression_summary_fmg_cube = tiledb.open(f"{corpus}/{EXPRESSION_SUMMARY_FMG_CUBE_NAME}")
        cell_counts_cube = tiledb.open(f"{corpus}/{CELL_COUNTS_CUBE_NAME}")
        dataset_to_gene_ids = json.load(open(f"{corpus}/{DATASET_TO_GENE_IDS_FILENAME}"))
    else:
        if corpus is None:
            corpus = load_snapshot()
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
        cell_type_orderings=None,
        primary_filter_dimensions=None,
    )
    q = WmgQuery(snapshot)
    query = q.expression_summary_fmg(criteria)
    cell_counts_query = q.cell_counts(criteria)

    # if group-by dimensions not provided, use the keys specified in the filter JSON
    if group_by_dims is None:
        # depluralize plural keys to match names in schema
        depluralized_keys = [i[:-1] if i[-1] == "s" else i for i in filters.keys()]
    else:
        depluralized_keys = [i[:-1] if i[-1] == "s" else i for i in group_by_dims]

    # group-by dimensions for expression summary must include genes
    gb_dims_es = ["gene_ontology_term_id"] + depluralized_keys
    # group-by dimensions for cell counts must include dataset ID
    # by convention, dataset_id will be first entry
    if "dataset_id" in depluralized_keys:
        depluralized_keys.remove("dataset_id")
    gb_dims = ["dataset_id"] + depluralized_keys

    # group-by and sum
    agg = query.groupby(gb_dims_es).sum()
    n_cells = cell_counts_query.groupby(gb_dims).sum()["n_total_cells"]

    keep_dataset_ids = group_by_dims is None and "dataset_id" in filters.keys()
    if genes is None:
        genes = list(agg.index.levels[0])

    t_n_cells_sum = _calculate_true_n_cells(n_cells, genes, dataset_to_gene_ids, keep_dataset_ids)
    return agg, t_n_cells_sum, genes


def _calculate_true_n_cells(n_cells, genes, dataset_to_gene_ids, keep_dataset_ids):
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

    keep_dataset_ids - boolean
        If True, then the dataset IDs are kept in the index of the returned
        population sizes. Otherwise, the dataset IDs are removed.
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
    if keep_dataset_ids:
        groups = list(n_cells.index)
    else:
        groups = list(zip(*[n_cells.index.get_level_values(i) for i in range(1, len(n_cells.index[0]))]))

    # sum up the cell count arrays across duplicate groups
    # (groups can be duplicate after excluding dataset_ids)
    t_n_cells_sum = {}
    for i, k in enumerate(groups):
        summer = t_n_cells_sum.get(k, np.zeros(t_n_cells.shape[1]))
        summer += t_n_cells[i]
        t_n_cells_sum[k] = summer
    return t_n_cells_sum


def _prepare_indices_and_metrics(target_filters, context_filters, corpus=None):
    """
    Compute all the necessary indices and metrics for the given target and context filters.

    Arguments
    ---------
    target_filters - dict
        Dictionary of filters for the target population

    context_filters - dict
        Dictionary of filters describing the context

    corpus - str or WmgSnapshot, optional, default None
        If string, it is the path to the snapshot.
        If WmgSnapshot, it is the snapshot object.
        If None, the snapshot will be fetched from AWS.

    Returns
    -------
    context_agg - DataFrame,
        aggregated values for the context query across the combinations of
        filters specified in target_filters

    target_agg - DataFrame,
        aggregated values for the target query across genes

    groups_context_uniq - list of tuples
        List of unique combinations of filter values (length N)

    n_target - np.ndarray
        1 x M array of cell counts for the target group
        (M = number of genes)

    n_context - np.ndarray
        N x M array of cell counts for the groups present in the context
        (N = number of groups, M = number of genes)

    target_data_nnz - np.ndarray
        1 x M array of number of nonzero values for each gene in the target group

    genes - list
        List of non-zero genes present in the target population

    groups_indices_context - list
        indices of the groups in `context_agg` corresponding to the index location of each group in
        `groups_context_uniq`.
        These are the row coordinates of the array in which the expression statistics will be filled
        for the context population.

    genes_indices_context - list
        indices of the genes in `context_agg` corresponding to the index location of each gene in
        `genes`.
        These are the column coordinates of the array in which the expression statistics will be filled
        for the context population.

    genes_indices_target - list
        indices of the genes in `target_agg` corresponding to the index location of each gene in
        `genes`.
        These are the column coordinates of the array in which the expression statistics will be filled
        for the target population.
    """
    context_agg, t_n_cells_sum_context, genes = _query_tiledb(
        context_filters, corpus=corpus, group_by_dims=list(target_filters.keys())
    )
    target_agg, t_n_cells_sum_target, _ = _query_tiledb(target_filters, corpus=corpus, genes=genes)
    groups_context_uniq = list(t_n_cells_sum_context.keys())

    target_agg = target_agg.groupby("gene_ontology_term_id").sum()

    genes_target = list(target_agg.index)

    genes_context = list(context_agg.index.get_level_values(0))
    groups_context = [i[1:] for i in context_agg.index]
    groups_indexer = pd.Series(index=groups_context_uniq, data=range(len(groups_context_uniq)))
    genes_indexer = pd.Series(index=genes, data=np.arange(len(genes)))

    groups_indices_context = list(groups_indexer[groups_context])
    genes_indices_context = list(genes_indexer[genes_context])
    genes_indices_target = list(genes_indexer[genes_target])

    target_data_nnz = np.zeros((1, len(genes)))
    target_data_nnz[0, genes_indices_target] = list(target_agg["nnz"])

    n_target = np.vstack(list(t_n_cells_sum_target.values()))
    n_context = np.vstack(list(t_n_cells_sum_context.values()))
    return (
        context_agg,
        target_agg,
        groups_context_uniq,
        n_target,
        n_context,
        target_data_nnz,
        genes,
        groups_indices_context,
        genes_indices_context,
        genes_indices_target,
    )


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


def _post_process_stats(genes, pvals, effects, nnz, test="ttest", min_num_expr_cells=25, percentile=0.15, n_markers=10):
    """
    Process and aggregate statistics into the output dictionary format.

    Arguments
    ---------
    genes - list
        List of genes corresponding to each p-value and effect size

    pvals - np.ndarray
        N x M array of p-values for each comparison and each gene

    effects - np.ndarray
        N x M array of effect sizes for each comparison and each gene

    nnz - np.ndarray
        1 x M array of number of nonzero expressions for target population (> 0)

    test - str, optional, default "ttest"
        The statistical test used ("ttest" or "binomtest").

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
    # aggregate
    effects = np.nanpercentile(effects, percentile * 100, axis=0)
    pvals = np.array([stats.combine_pvalues(x[np.invert(np.isnan(x))] + 1e-300)[-1] for x in pvals.T])
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


def _get_markers_ttest(target_filters, context_filters, corpus=None, n_markers=10, percentile=0.15):
    """
    Calculate marker genes using the t-test.

    Arguments
    ---------
    target_filters - dict
        Dictionary of filters for the target population

    context_filters - dict
        Dictionary of filters describing the context

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
    (
        context_agg,
        target_agg,
        groups_context_uniq,
        n_target,
        n_context,
        target_data_nnz,
        genes,
        groups_indices_context,
        genes_indices_context,
        genes_indices_target,
    ) = _prepare_indices_and_metrics(target_filters, context_filters, corpus=corpus)

    target_data_sum = np.zeros((1, len(genes)))
    target_data_sum[0, genes_indices_target] = list(target_agg["sum"])
    target_data_sumsq = np.zeros((1, len(genes)))
    target_data_sumsq[0, genes_indices_target] = list(target_agg["sqsum"])

    context_data_sum = np.zeros((len(groups_context_uniq), len(genes)))
    context_data_sum[groups_indices_context, genes_indices_context] = list(context_agg["sum"])
    context_data_sumsq = np.zeros((len(groups_context_uniq), len(genes)))
    context_data_sumsq[groups_indices_context, genes_indices_context] = list(context_agg["sqsum"])

    pvals, effects = _run_ttest(
        target_data_sum, target_data_sumsq, n_target, context_data_sum, context_data_sumsq, n_context
    )
    return _post_process_stats(
        genes,
        pvals,
        effects,
        target_data_nnz,
        test="ttest",
        n_markers=n_markers,
        percentile=percentile,
    )


def _run_binom(nnz_thr1, n1, nnz_thr2, n2):
    """
    Calculate binomial test statistics.

    Arguments
    ---------
    nnz_thr1 - np.ndarray (1 x M)
        Array of number of nonzero expressions greater than a threshold (1) for target pop for each gene
    n1 - np.ndarray (1 x M)
        Number of cells in target population for each gene
    nnz_thr2 - np.ndarray (N x M)
        Array of number of nonzero expressions greater than a threshold (1) for context pops for each gene
    n2 - np.ndarray (N x M)
        Number of cells in context populations for each gene

    Returns
    -------
    pvals_adj - np.ndarray
        adjusted p-values for each comparison for each gene
    effects - np.ndarray
        effect sizes for each comparison for each gene
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        p_query = n1 / (n1 + n2)
        size = nnz_thr1 + nnz_thr2
        pvals = stats.binom.sf(nnz_thr1 + 1, size + 2, p_query)
        mean_n = (n1 + n2) / 2
        pn1 = n1 / mean_n
        pn2 = n2 / mean_n
        effects = np.log2((nnz_thr1 + pn1) / (n1 + 2 * pn1)) - np.log2((nnz_thr2 + pn2) / (n2 + 2 * pn2))

    return pvals, effects


def _get_markers_binomtest(target_filters, context_filters, corpus=None, n_markers=10, percentile=0.8):
    """
    Calculate marker genes using the binomial test.

    Arguments
    ---------
    target_filters - dict
        Dictionary of filters for the target population

    context_filters - dict
        Dictionary of filters describing the context

    corpus - str or WmgSnapshot, optional, default None
        If string, it is the path to the snapshot.
        If WmgSnapshot, it is the snapshot object.
        If None, the snapshot will be fetched from AWS.

    n_markers - int, optional, default 10
        Number of top markers to return. If None, all marker genes with effect size > 0 are returned.

    percentile - float, optional, default 0.8
        The percentile of effect sizes to select as the representative effect size.

    Returns
    -------
    Dictionary of marker genes as keys and a dictionary of p-value and effect size as values.
    """

    (
        context_agg,
        target_agg,
        groups_context_uniq,
        n_target,
        n_context,
        target_data_nnz,
        genes,
        groups_indices_context,
        genes_indices_context,
        genes_indices_target,
    ) = _prepare_indices_and_metrics(target_filters, context_filters, corpus=corpus)
    target_data_nnz_thr = np.zeros((1, len(genes)))
    target_data_nnz_thr[0, genes_indices_target] = list(target_agg["nnz_thr"])
    context_data_nnz_thr = np.zeros((len(groups_context_uniq), len(genes)))
    context_data_nnz_thr[groups_indices_context, genes_indices_context] = list(context_agg["nnz_thr"])

    pvals, effects = _run_binom(target_data_nnz_thr, n_target, context_data_nnz_thr, n_context)

    return _post_process_stats(
        genes,
        pvals,
        effects,
        target_data_nnz,
        test="binomtest",
        n_markers=n_markers,
        percentile=percentile,
    )


def get_markers(target_filters, context_filters, corpus=None, test="ttest", n_markers=10, percentile=0.15):
    """
    Calculate marker genes using the t-test.

    Arguments
    ---------
    target_filters - dict
        Dictionary of filters for the target population

    context_filters - dict
        Dictionary of filters describing the context

    corpus - str or WmgSnapshot, optional, default None
        If string, it is the path to the snapshot.
        If WmgSnapshot, it is the snapshot object.
        If None, the snapshot will be fetched from AWS.

    test - str, optional, default "ttest"
        The statistical test to be used ("ttest" or "binomtest").

    n_markers - int, optional, default 10
        Number of top markers to return. If None, all marker genes with effect size > 0 are returned.

    percentile - float, optional, default 0.15
        The percentile of effect sizes to select as the representative effect size.

    Returns
    -------
    Dictionary of marker genes as keys and a dictionary of p-value and effect size as values.
    """

    if test == "ttest":
        return _get_markers_ttest(
            target_filters,
            context_filters,
            corpus=corpus,
            n_markers=n_markers,
            percentile=percentile,
        )
    elif test == "binomtest":
        return _get_markers_binomtest(
            target_filters,
            context_filters,
            corpus=corpus,
            n_markers=n_markers,
            percentile=percentile,
        )
    else:
        raise ValueError(f"Test {test} not supported.")
