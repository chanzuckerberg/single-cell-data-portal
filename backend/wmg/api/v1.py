from collections import defaultdict
from typing import Dict, List, Any, Iterable, Tuple
from math import isnan

import connexion
from flask import jsonify
from pandas import DataFrame

from backend.common.entities import Dataset
from backend.common.utils.db_session import db_session_manager
from backend.wmg.data.ontology_labels import ontology_term_label, gene_term_label
from backend.wmg.api.query import (
    WmgQuery,
    WmgQueryCriteria,
)
from backend.wmg.data.snapshot import load_snapshot, WmgSnapshot

# TODO: add cache directives: no-cache (i.e. revalidate); impl etag
#  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data
#  -portal/2132


def primary_filter_dimensions():
    snapshot: WmgSnapshot = load_snapshot()
    return jsonify(snapshot.primary_filter_dimensions)


def query():
    request = connexion.request.json
    criteria = WmgQueryCriteria(**request["filter"])

    snapshot: WmgSnapshot = load_snapshot()
    query = WmgQuery(snapshot)
    expression_summary = query.expression_summary(criteria)
    cell_counts = query.cell_counts(criteria)
    dot_plot_matrix_df, cell_counts_cell_type_agg = get_dot_plot_data(expression_summary, cell_counts)

    include_filter_dims = request.get("include_filter_dims", False)

    response_filter_dims_values = (
        build_filter_dims_values(criteria, query, expression_summary) if include_filter_dims else {}
    )
    return jsonify(
        dict(
            snapshot_id=snapshot.snapshot_identifier,
            expression_summary=build_expression_summary(dot_plot_matrix_df),
            term_id_labels=dict(
                genes=build_gene_id_label_mapping(criteria.gene_ontology_term_ids),
                cell_types=build_ordered_cell_types_by_tissue(
                    cell_counts, cell_counts_cell_type_agg.T, snapshot.cell_type_orderings
                ),
            ),
            filter_dims=response_filter_dims_values,
        )
    )


def fmg_query():
    request = connexion.request.json
    target_filters = FmgQueryCriteria(**request["target_filters"])
    context_filters = FmgQueryCriteria(**request["context_filters"])

    snapshot: WmgSnapshot = load_snapshot()
    query = WmgQuery(snapshot)
    
    target_query = query.expression_summary_fmg(target_filters)
    context_query = query.expression_summary_fmg(context_filters)
    target_cell_counts_query = query.cell_counts(target_filters)
    context_cell_counts_query = query.cell_counts(context_filters)

    gb_dims_es = ['gene_ontology_term_id']+list(target_filters.keys())
    context_agg = context_query.groupby(gb_dims_es).sum()
    target_agg = target_query.groupby(gb_dims_es).sum()

    genes = list(target_agg.index.levels[0])
    genes_indexer = pd.Series(index=genes,data=range(len(genes)))

    filt = np.in1d(np.array(list(context_agg.index.get_level_values(0))),genes)
    context_agg = context_agg[filt]

    gb_dims = ['dataset_id']+list(target_filters.keys())
    n_cells_target = target_cell_counts_query.groupby(gb_dims).sum()['n_cells']
    n_cells_context = context_cell_counts_query.groupby(gb_dims).sum()['n_cells']

    dataset_to_genes_present = {}
    for i in set(n_cells_target.index.levels[0]).union(n_cells_context.index.levels[0]):
        dataset_to_genes_present[i] = list(set(dataset_to_genes[i]).intersection(genes))
        
    desired_levels = [i.name for i in target_agg.index.levels][1:]

    def fill_in_n_cells(n_cells, agg):
        t_n_cells = np.zeros((n_cells.shape[0],len(genes)))
        
        for i, index in enumerate(n_cells.index):
            n = n_cells[index]
            dataset_id = index[0]
            t_n_cells[i,genes_indexer[dataset_to_genes_present[dataset_id]]]=n
        
        ixs = []
        level_names = [i.name for i in n_cells.index.levels]
        for d in desired_levels:
            for i,level in enumerate(level_names):
                if level == d:
                    ixs.append(i)
        
        groups = list(zip(*[n_cells.index.get_level_values(i) for i in ixs]))

        t_n_cells_sum = {}
        for i,k in enumerate(groups):
            summer = t_n_cells_sum.get(k, np.zeros(t_n_cells.shape[1]))
            summer += t_n_cells[i]
            t_n_cells_sum[k] = summer
            
        ns = []
        for index in agg.index:
            try:
                ns.append(
                    t_n_cells_sum[tuple(index[1:])][genes_indexer[index[0]]]
                )
            except KeyError:
                ns.append(0)
        agg['n_cells'] = ns
        return list(t_n_cells_sum.keys()), t_n_cells_sum
        
    groups_target_uniq, t_n_cells_sum_target = fill_in_n_cells(n_cells_target, target_agg)
    groups_context_uniq, t_n_cells_sum_context = fill_in_n_cells(n_cells_context, context_agg)

    target_agg = target_agg.groupby("gene_ontology_term_id").sum()

    for df in [target_agg,context_agg]:
        df['mean'] = df['sum'] / df['n_cells']
        df['mean_sq'] = df['sqsum'] / df['n_cells']

    genes_target = genes
    genes_uniq = genes

    genes_context = list(context_agg.index.get_level_values(0))
    groups_context = [i[1:] for i in context_agg.index]

    nc_context = np.array(context_agg['n_cells'])

    groups_indexer = pd.Series(index=groups_context_uniq,data=range(len(groups_context_uniq)))
    genes_indexer = pd.Series(index=genes_uniq,data=np.arange(len(genes_uniq)))

    groups_indices_context = list(groups_indexer[groups_context])
    genes_indices_context = list(genes_indexer[genes_context])
    genes_indices_target = list(genes_indexer[genes_target])

    #context_data_nnz = np.zeros((len(x_context_uniq),len(y_uniq)))
    #context_data_nnz[x_indices_context,y_indices_context]=list(context_agg['nnz'])

    context_data_mean = np.zeros((len(groups_context_uniq),len(genes_uniq)))
    context_data_mean[groups_indices_context,genes_indices_context]=list(context_agg['mean'])

    context_data_meansq = np.zeros((len(groups_context_uniq),len(genes_uniq)))
    context_data_meansq[groups_indices_context,genes_indices_context]=list(context_agg['mean_sq'])
    context_data_var = context_data_meansq - context_data_mean**2
    context_data_var[context_data_var<0]=0

    target_data_nnz = np.zeros((1,len(genes_uniq)))
    target_data_nnz[0,genes_indices_target]=list(target_agg['nnz'])

    target_data_mean = np.zeros((1,len(genes_uniq)))
    target_data_mean[0,genes_indices_target]=list(target_agg['mean'])

    target_data_meansq = np.zeros((1,len(genes_uniq)))
    target_data_meansq[0,genes_indices_target]=list(target_agg['mean_sq'])
    target_data_var = target_data_meansq - target_data_mean**2
    target_data_var[target_data_var<0]=0

    n1 = np.vstack(list(t_n_cells_sum_target.values()))
    n2 = np.vstack(list(t_n_cells_sum_context.values()))

    with np.errstate(divide="ignore", invalid="ignore"):
        target_data_var_n = target_data_var / n1
        context_data_var_n = context_data_var / n2
        sum_data_var_n = target_data_var_n + context_data_var_n    
        dof = sum_data_var_n**2 / (target_data_var_n**2 / (n1 - 1) + context_data_var_n**2 / (n2 - 1))
        dof[np.isnan(dof)] = 1
        tscores = (target_data_mean - context_data_mean) / np.sqrt(sum_data_var_n)
        effects = (target_data_mean - context_data_mean) / np.sqrt( ((n1 - 1) * target_data_var + (n2 - 1) * context_data_var) / (n1+n2-1))
        
    tscores[np.isnan(tscores)] = 0
    effects[np.isnan(effects)] = 0
    pvals = stats.t.sf(tscores, dof)
    zero_out = target_data_nnz.flatten() < 50

    tscores[:,zero_out]=0
    effects[:,zero_out]=0
    pvals[:,zero_out]=1

    pvals_adj = pvals * target_data_var.size
    pvals_adj[pvals_adj > 1] = 1  # cap adjusted p-value at 1

    # aggregate
    tscores = np.sort(tscores,axis=0)[:N_BOTTOM_COMPARISONS].mean(0)
    effects = np.sort(effects,axis=0)[:N_BOTTOM_COMPARISONS].mean(0)
    #todo: fix p-value aggregation
    pvals_adj = np.sort(pvals_adj,axis=0)[:N_BOTTOM_COMPARISONS].mean(0)

    markers = np.array(genes_uniq)[np.argsort(-effects)[:N_MARKERS]]
    p = pvals_adj[np.argsort(-effects)[:N_MARKERS]]    
    effects = effects[np.argsort(-effects)][:N_MARKERS]
    
    return jsonify(dict(zip(markers,list(zip(p,effects)))))

# TODO: Read this from generated data artifact instead of DB.
#  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data
#  -portal/2086. This code is without a unit test, but we are intending to replace it.
def fetch_datasets_metadata(dataset_ids: Iterable[str]) -> List[Dict]:
    # We return a DTO because the db entity can't access its attributes after the db session ends,
    # and we want to keep session management out of the calling method

    with db_session_manager() as session:

        def get_dataset(dataset_id_):
            dataset = Dataset.get(session, dataset_id_)
            if dataset is None:
                # Handle possible missing dataset due to db state evolving past wmg snapshot
                return dict(id=dataset_id_, label="", collection_id="", collection_label="")
            return dict(
                id=dataset.id,
                label=dataset.name,
                collection_id=dataset.collection.id,
                collection_label=dataset.collection.name,
            )

        return [get_dataset(dataset_id) for dataset_id in dataset_ids]


def find_dim_option_values(criteria: Dict, query: WmgQuery, dimension: str) -> set:
    """Find values for the specified dimension that satisfy the given filtering criteria,
    ignoring any criteria specified for the given dimension."""
    filter_options_criteria = criteria.copy(update={dimension + "s": []}, deep=True)
    # todo can we query cell_counts for a performance gain?
    query_result = query.expression_summary(filter_options_criteria)
    filter_dims = query_result.groupby(dimension).groups.keys()
    return filter_dims


def build_filter_dims_values(criteria: WmgQueryCriteria, query: WmgQuery, expression_summary: DataFrame) -> Dict:
    dims = {
        "dataset_id": "",
        "disease_ontology_term_id": "",
        "sex_ontology_term_id": "",
        "development_stage_ontology_term_id": "",
        "ethnicity_ontology_term_id": "",
    }
    for dim in dims:
        if len(criteria.dict()[dim + "s"]) == 0:
            dims[dim] = expression_summary.groupby(dim).groups.keys()
        else:
            dims[dim] = find_dim_option_values(criteria, query, dim)

    response_filter_dims_values = dict(
        datasets=fetch_datasets_metadata(dims["dataset_id"]),
        disease_terms=build_ontology_term_id_label_mapping(dims["disease_ontology_term_id"]),
        sex_terms=build_ontology_term_id_label_mapping(dims["sex_ontology_term_id"]),
        development_stage_terms=build_ontology_term_id_label_mapping(dims["development_stage_ontology_term_id"]),
        ethnicity_terms=build_ontology_term_id_label_mapping(dims["ethnicity_ontology_term_id"]),
    )

    return response_filter_dims_values


def build_expression_summary(query_result: DataFrame) -> dict:
    # Create nested dicts with gene_ontology_term_id, tissue_ontology_term_id keys, respectively
    structured_result: Dict[str, Dict[str, List[Dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
    for row in query_result.itertuples(index=False):
        structured_result[row.gene_ontology_term_id][row.tissue_ontology_term_id].append(
            dict(
                id=row.cell_type_ontology_term_id,
                n=row.nnz,
                me=row.sum / row.nnz,
                pc=row.nnz / row.n_cells_cell_type,
                tpc=row.nnz / row.n_cells_tissue,
            )
        )
    return structured_result


def agg_cell_type_counts(cell_counts: DataFrame) -> DataFrame:
    # Aggregate cube data by tissue, cell type
    cell_counts_cell_type_agg = cell_counts.groupby(
        ["tissue_ontology_term_id", "cell_type_ontology_term_id"], as_index=True
    ).sum()
    cell_counts_cell_type_agg.rename(columns={"n_total_cells": "n_cells_cell_type"}, inplace=True)
    return cell_counts_cell_type_agg


def agg_tissue_counts(cell_counts: DataFrame) -> DataFrame:
    # Aggregate cube data by tissue
    cell_counts_tissue_agg = cell_counts.groupby(["tissue_ontology_term_id"], as_index=True).sum()
    cell_counts_tissue_agg.rename(columns={"n_total_cells": "n_cells_tissue"}, inplace=True)
    return cell_counts_tissue_agg


def get_dot_plot_data(query_result: DataFrame, cell_counts: DataFrame) -> Tuple[DataFrame, DataFrame]:
    # Get the dot plot matrix dataframe and aggregated cell counts per cell type
    cell_counts_cell_type_agg = agg_cell_type_counts(cell_counts)
    cell_counts_tissue_agg = agg_tissue_counts(cell_counts)
    dot_plot_matrix_df = build_dot_plot_matrix(query_result, cell_counts_cell_type_agg, cell_counts_tissue_agg)
    return dot_plot_matrix_df, cell_counts_cell_type_agg


def build_dot_plot_matrix(
    query_result: DataFrame, cell_counts_cell_type_agg: DataFrame, cell_counts_tissue_agg: DataFrame
) -> DataFrame:
    # Aggregate cube data by gene, tissue, cell type
    expr_summary_agg = query_result.groupby(
        ["gene_ontology_term_id", "tissue_ontology_term_id", "cell_type_ontology_term_id"], as_index=False
    ).sum()
    return expr_summary_agg.join(
        cell_counts_cell_type_agg, on=["tissue_ontology_term_id", "cell_type_ontology_term_id"], how="left"
    ).join(cell_counts_tissue_agg, on=["tissue_ontology_term_id"], how="left")


def build_gene_id_label_mapping(gene_ontology_term_ids: List[str]) -> List[dict]:
    return [
        {gene_ontology_term_id: gene_term_label(gene_ontology_term_id)}
        for gene_ontology_term_id in gene_ontology_term_ids
    ]


def build_ontology_term_id_label_mapping(ontology_term_ids: Iterable[str]) -> List[dict]:
    return [{ontology_term_id: ontology_term_label(ontology_term_id)} for ontology_term_id in ontology_term_ids]


def build_ordered_cell_types_by_tissue(
    cell_counts: DataFrame,
    cell_counts_cell_type_agg_T: DataFrame,
    cell_type_orderings: DataFrame,
) -> Dict[str, List[Dict[str, str]]]:
    distinct_tissues_cell_types: DataFrame = cell_counts.groupby(
        ["tissue_ontology_term_id", "cell_type_ontology_term_id"], as_index=False
    ).first()[["tissue_ontology_term_id", "cell_type_ontology_term_id", "n_total_cells"]]

    joined = cell_type_orderings.merge(
        distinct_tissues_cell_types, on=["tissue_ontology_term_id", "cell_type_ontology_term_id"], how="left"
    )

    # Updates depths based on the rows that need to be removed
    joined = build_ordered_cell_types_by_tissue_update_depths(joined)

    # Remove cell types without counts
    joined = joined[joined["n_total_cells"].notnull()]

    structured_result: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in joined.itertuples(index=False):
        structured_result[row.tissue_ontology_term_id].append(
            {
                "cell_type_ontology_term_id": row.cell_type_ontology_term_id,
                "cell_type": ontology_term_label(row.cell_type_ontology_term_id),
                "total_count": int(
                    cell_counts_cell_type_agg_T[row.tissue_ontology_term_id][row.cell_type_ontology_term_id][
                        "n_cells_cell_type"
                    ]
                ),
                "depth": row.depth,
            }
        )

    return structured_result


def build_ordered_cell_types_by_tissue_update_depths(x: DataFrame):
    """
    Updates the depths of the cell ontology tree based on cell types that have to be removed
    because they have 0 counts
    """

    depth_col = x.columns.get_loc("depth")
    n_cells_col = x.columns.get_loc("n_total_cells")

    x["depth"] = x["depth"].astype("int")

    for i in range(len(x)):
        if isnan(x.iloc[i, n_cells_col]):
            original_depth = x.iloc[i, depth_col]
            for j in range(i + 1, len(x)):
                if original_depth < x.iloc[j, depth_col]:
                    x.iloc[j, depth_col] -= 1
                else:
                    break

    return x
