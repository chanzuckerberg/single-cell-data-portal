import pandas as pd
import numpy as np
from scipy import stats
import json
from functools import lru_cache, wraps
from backend.wmg.api.query import WmgQuery, FmgQueryCriteria
from backend.wmg.data.snapshot import load_snapshot, WmgSnapshot

def _make_hashable(func):
    class HDict(dict):
        def __hash__(self):
            return hash(json.dumps(self))

    class HList(list):
        def __hash__(self):
            return hash(tuple(self))

    @wraps(func)
    def wrapped(*args, **kwargs):
        args = tuple([HDict(arg) if isinstance(arg, dict) else arg for arg in args])
        kwargs = {k: HDict(v) if isinstance(v, dict) else v for k, v in kwargs.items()}

        args = tuple([HList(arg) if isinstance(arg, list) else arg for arg in args])
        kwargs = {k: HList(v) if isinstance(v, list) else v for k, v in kwargs.items()}
        return func(*args, **kwargs)

    return wrapped


@_make_hashable
@lru_cache(maxsize=None)
def _query_tiledb(filters, group_by_dims=None, genes=None):
    criteria = FmgQueryCriteria(**filters)
    snapshot: WmgSnapshot = load_snapshot()
    Q = WmgQuery(snapshot)
    query = Q.expression_summary_fmg(criteria)
    cell_counts_query = Q.cell_counts(criteria)

    if group_by_dims is None:
        gb_dims_es = ["gene_ontology_term_id"] + list(filters.keys())
        gb_dims = ["dataset_id"] + list(filters.keys())
    else:
        gb_dims_es = ["gene_ontology_term_id"] + group_by_dims
        gb_dims = ["dataset_id"] + group_by_dims

    agg = query.groupby(gb_dims_es).sum()
    n_cells = cell_counts_query.groupby(gb_dims).sum()["n_cells"]

    if group_by_dims is None:
        desired_levels = [i.name for i in agg.index.levels][1:]
    else:
        desired_levels = group_by_dims

    if genes is None:
        genes = list(agg.index.levels[0])

    genes_indexer = pd.Series(index=genes, data=range(len(genes)))

    t_n_cells = np.zeros((n_cells.shape[0], len(genes)))

    for i, index in enumerate(n_cells.index):
        n = n_cells[index]
        dataset_id = index[0]
        present = list(set(snapshot.dataset_to_gene_ids[dataset_id]).intersection(genes))
        t_n_cells[i, genes_indexer[present]] = n

    ixs = []
    level_names = [i.name for i in n_cells.index.levels]
    for d in desired_levels:
        for i, level in enumerate(level_names):
            if level == d:
                ixs.append(i)

    groups = list(zip(*[n_cells.index.get_level_values(i) for i in ixs]))

    t_n_cells_sum = {}
    for i, k in enumerate(groups):
        summer = t_n_cells_sum.get(k, np.zeros(t_n_cells.shape[1]))
        summer += t_n_cells[i]
        t_n_cells_sum[k] = summer

    return agg, t_n_cells_sum, genes


@_make_hashable
def _prepare_indices_and_metrics(target_filters, context_filters):
    context_agg, groups_context_uniq, t_n_cells_sum_context, genes = _query_tiledb(
        context_filters, group_by_dims=list(target_filters.keys())
    )
    target_agg, t_n_cells_sum_target, _ = _query_tiledb(target_filters, genes=genes)

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


def _post_process_stats(genes, pvals, effects, nnz, min_num_expr_cells=25, p_bottom_comparisons=0.1, n_markers=10):
    zero_out = nnz.flatten() < min_num_expr_cells
    effects[:, zero_out] = 0
    pvals[:, zero_out] = 0
    # aggregate
    n_bottom_comparisons = int(p_bottom_comparisons * effects.shape[0]) + 1
    effects = np.sort(effects, axis=0)[:n_bottom_comparisons].mean(0)
    # todo: fix p-value aggregation
    pvals = np.sort(pvals, axis=0)[:n_bottom_comparisons].mean(0)
    markers = np.array(genes)[np.argsort(-effects)[:n_markers]]
    p = pvals[np.argsort(-effects)[:n_markers]]
    effects = effects[np.argsort(-effects)[:n_markers]]
    statistics = []
    for i in range(len(p)):
        pi = p[i]
        ei = effects[i]
        statistics.append({"p_value": pi, "effect_size": ei})
    return dict(zip(list(markers), statistics))


def _get_markers_ttest(target_filters, context_filters, n_markers=10, p_bottom_comparisons=0.1):
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
    ) = _prepare_indices_and_metrics(target_filters, context_filters)

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
        genes, pvals, effects, target_data_nnz, n_markers=n_markers, p_bottom_comparisons=p_bottom_comparisons
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


def _get_markers_binomtest(target_filters, context_filters, n_markers=10, p_bottom_comparisons=0.8):
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
    ) = _prepare_indices_and_metrics(target_filters, context_filters)
    target_data_nnz_thr = np.zeros((1, len(genes)))
    target_data_nnz_thr[0, genes_indices_target] = list(target_agg["nnz_thr"])
    context_data_nnz_thr = np.zeros((len(groups_context_uniq), len(genes)))
    context_data_nnz_thr[groups_indices_context, genes_indices_context] = list(context_agg["nnz_thr"])

    pvals, effects = _run_binom(target_data_nnz_thr, n_target, context_data_nnz_thr, n_context)

    return _post_process_stats(
        genes, pvals, effects, target_data_nnz, n_markers=n_markers, p_bottom_comparisons=p_bottom_comparisons
    )


def get_markers(target_filters, context_filters, test="ttest", n_markers=10, p_bottom_comparisons=0.1):
    if test == "ttest":
        return _get_markers_ttest(
            target_filters, context_filters, n_markers=n_markers, p_bottom_comparisons=p_bottom_comparisons
        )
    elif test == "binom":
        return _get_markers_binomtest(
            target_filters, context_filters, n_markers=n_markers, p_bottom_comparisons=p_bottom_comparisons
        )
    else:
        raise ValueError(f"Test {test} not supported.")
