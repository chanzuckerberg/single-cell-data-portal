import pandas as pd
import numpy as np
from scipy import stats
import json
from functools import lru_cache, wraps
from backend.wmg.api.query import WmgQuery, FmgQueryCriteria
from backend.wmg.data.snapshot import load_snapshot, WmgSnapshot


def _make_hashable(func):
    """
    Implicitly convert unhashable data structures (list and dict) to hashable strings
    for memoization.
    """
    class HDict(dict):
        def __hash__(self):
            return hash(json.dumps(self))

    class HList(list):
        def __hash__(self):
            return hash(tuple(self))

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

@_make_hashable
@lru_cache(maxsize=None)
def _query_tiledb(filters, corpus_path=None, group_by_dims=None, genes=None):
    """
    Query tileDB and return the following results:
    - agg: query result aggregated by filter keys
    - t_n_cells_sum: a dictionary of 1D numpy arrays with length # of genes
        each value is the true number of cells present for the combination of 
        filters (keys).
        Ex: {(cell_type1,tissue1,organism1): np.array([0,100,100,0,80])}
    - genes: a list of all genes with nonzero expression found in the query result
    """

    criteria = FmgQueryCriteria(**filters)
    snapshot: WmgSnapshot = load_snapshot(corpus_path=corpus_path)
    Q = WmgQuery(snapshot)
    query = Q.expression_summary_fmg(criteria)
    cell_counts_query = Q.cell_counts(criteria)

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
    if 'dataset_id' in depluralized_keys:
        depluralized_keys.remove('dataset_id')
    gb_dims = ["dataset_id"] + depluralized_keys

    # group-by and sum
    agg = query.groupby(gb_dims_es).sum()
    n_cells = cell_counts_query.groupby(gb_dims).sum()["n_total_cells"]

    if group_by_dims is None:
        desired_levels = [i.name for i in agg.index.levels][1:]
    else:
        desired_levels = depluralized_keys

    if genes is None:
        genes = list(agg.index.levels[0])

    genes_indexer = pd.Series(index=genes, data=range(len(genes)))

    # fill in an array stratified across desired filter combinations + dataset ID
    # with the number of cells in each group for genes that are PRESENT in the dataset
    t_n_cells = np.zeros((n_cells.shape[0], len(genes)))
    for i, index in enumerate(n_cells.index):
        n = n_cells[index]
        dataset_id = index[0]
        present = list(set(snapshot.dataset_to_gene_ids[dataset_id]).intersection(genes))
        t_n_cells[i, genes_indexer[present]] = n

    # get the tuples of filter values in the correct order defined in `desired_levels`
    # this is to exclude dataset_ids from the groups
    ixs = []
    level_names = [i.name for i in n_cells.index.levels]
    for d in desired_levels:
        for i, level in enumerate(level_names):
            if level == d:
                ixs.append(i)
    groups = list(zip(*[n_cells.index.get_level_values(i) for i in ixs]))

    # sum up the cell count arrays across duplicate groups 
    # (groups can be duplicate after excluding dataset_ids)
    t_n_cells_sum = {}
    for i, k in enumerate(groups):
        summer = t_n_cells_sum.get(k, np.zeros(t_n_cells.shape[1]))
        summer += t_n_cells[i]
        t_n_cells_sum[k] = summer

    return agg, t_n_cells_sum, genes


def _prepare_indices_and_metrics(target_filters, context_filters, corpus_path=None):
    context_agg, t_n_cells_sum_context, genes = _query_tiledb(
        context_filters, corpus_path=corpus_path, group_by_dims=list(target_filters.keys())
    )
    target_agg, t_n_cells_sum_target, _ = _query_tiledb(target_filters, corpus_path=corpus_path, genes=genes)
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
        groups_context_uniq,
        target_agg,
        n_target,
        n_context,
        target_data_nnz,
        genes,
        groups_indices_context,
        genes_indices_context,
        genes_indices_target,
    )


def _run_ttest(sum1, sumsq1, n1, sum2, sumsq2, n2):
    with np.errstate(divide="ignore", invalid="ignore"):
        mean1 = sum1 / n1
        meansq1 = sumsq1 / n1

        mean2 = sum2 / n2
        meansq2 = sumsq2 / n2

        for i in [mean1, meansq1, mean2, meansq2]:
            i[np.isnan(i)] = 0

        var1 = meansq1 - mean1**2
        var1[var1 < 0] = 0
        var2 = meansq2 - mean2**2
        var2[var2 < 0] = 0

        var1_n = var1 / n1
        var2_n = var2 / n2
        sum_var_n = var1_n + var2_n
        dof = sum_var_n**2 / (var1_n**2 / (n1 - 1) + var2_n**2 / (n2 - 1))
        dof[np.isnan(dof)] = 1
        tscores = (mean1 - mean2) / np.sqrt(sum_var_n)
        effects = (mean1 - mean2) / np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 1))

    tscores[np.isnan(tscores)] = 0
    effects[np.isnan(effects)] = 0
    pvals = stats.t.sf(tscores, dof)
    pvals_adj = pvals * sum1.size
    pvals_adj[pvals_adj > 1] = 1  # cap adjusted p-value at 1
    return pvals_adj, effects


def _post_process_stats(
    genes, pvals, effects, nnz, test="ttest", min_num_expr_cells=25, p_bottom_comparisons=0.1, n_markers=10
):
    zero_out = nnz.flatten() < min_num_expr_cells
    effects[:, zero_out] = 0
    pvals[:, zero_out] = 0
    # aggregate
    n_bottom_comparisons = int(p_bottom_comparisons * effects.shape[0]) + 1
    effects = np.sort(effects, axis=0)[:n_bottom_comparisons].mean(0)
    # todo: fix p-value aggregation
    pvals = np.sort(pvals, axis=0)[:n_bottom_comparisons].mean(0)
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
        if ei > 0:
            statistics.append({f"p_value_{test}": pi, f"effect_size_{test}": ei})
            final_markers.append(markers[i])
    return dict(zip(list(final_markers), statistics))


def _get_markers_ttest(target_filters, context_filters, corpus_path=None, n_markers=10, p_bottom_comparisons=0.1):
    (
        context_agg,
        groups_context_uniq,
        target_agg,
        n_target,
        n_context,
        target_data_nnz,
        genes,
        groups_indices_context,
        genes_indices_context,
        genes_indices_target,
    ) = _prepare_indices_and_metrics(target_filters, context_filters, corpus_path=corpus_path)

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
        p_bottom_comparisons=p_bottom_comparisons,
    )


def _run_binom(nnz_thr1, n1, nnz_thr2, n2):
    with np.errstate(divide="ignore", invalid="ignore"):
        p_query = n1 / (n1 + n2)
        size = nnz_thr1 + nnz_thr2
        pvals = stats.binom.sf(nnz_thr1, size, p_query)
        mean_n = (n1 + n2) / 2
        pn1 = n1 / mean_n
        pn2 = n2 / mean_n
        effects = np.log2((nnz_thr1 + pn1) / (n1 + 2 * pn1)) - np.log2((nnz_thr2 + pn2) / (n2 + 2 * pn2))

    effects[np.isnan(effects)] = 0

    pvals_adj = pvals * nnz_thr1.size
    pvals_adj[pvals_adj > 1] = 1  # cap adjusted p-value at 1
    return pvals_adj, effects


def _get_markers_binomtest(target_filters, context_filters, corpus_path=None, n_markers=10, p_bottom_comparisons=0.8):
    (
        context_agg,
        groups_context_uniq,
        target_agg,
        n_target,
        n_context,
        target_data_nnz,
        genes,
        groups_indices_context,
        genes_indices_context,
        genes_indices_target,
    ) = _prepare_indices_and_metrics(target_filters, context_filters, corpus_path=corpus_path)
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
        p_bottom_comparisons=p_bottom_comparisons,
    )


def get_markers(
    target_filters, context_filters, corpus_path=None, test="ttest", n_markers=10, p_bottom_comparisons=0.1
):
    if test == "ttest":
        return _get_markers_ttest(
            target_filters,
            context_filters,
            corpus_path=corpus_path,
            n_markers=n_markers,
            p_bottom_comparisons=p_bottom_comparisons,
        )
    elif test == "binomtest":
        return _get_markers_binomtest(
            target_filters,
            context_filters,
            corpus_path=corpus_path,
            n_markers=n_markers,
            p_bottom_comparisons=p_bottom_comparisons,
        )
    else:
        raise ValueError(f"Test {test} not supported.")
