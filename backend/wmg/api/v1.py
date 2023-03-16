from collections import defaultdict
from typing import Any, Dict, Iterable, List, Tuple

import connexion
from flask import jsonify
from pandas import DataFrame
from server_timing import Timing as ServerTiming

from backend.wmg.data.ontology_labels import gene_term_label, ontology_term_label
from backend.wmg.data.query import (
    MarkerGeneQueryCriteria,
    WmgFiltersQueryCriteria,
    WmgQuery,
    WmgQueryCriteria,
    retrieve_top_n_markers,
)
from backend.wmg.data.rollup import rollup_across_cell_type_descendants
from backend.wmg.data.schemas.cube_schema import expression_summary_non_indexed_dims
from backend.wmg.data.snapshot import WmgSnapshot, load_snapshot
from backend.wmg.data.utils import depluralize, find_all_dim_option_values, find_dim_option_values

# TODO: add cache directives: no-cache (i.e. revalidate); impl etag
#  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data
#  -portal/2132

DEFAULT_GROUP_BY_TERMS = ["tissue_ontology_term_id", "cell_type_ontology_term_id"]


def primary_filter_dimensions():
    snapshot: WmgSnapshot = load_snapshot()
    return jsonify(snapshot.primary_filter_dimensions)


def query():
    request = connexion.request.json
    is_rollup = request.get("is_rollup", True)
    compare = request.get("compare", None)

    if compare:
        compare = find_dimension_id_from_compare(compare)

    criteria = WmgQueryCriteria(**request["filter"])

    with ServerTiming.time("query and build response"):
        snapshot: WmgSnapshot = load_snapshot()
        q = WmgQuery(snapshot)
        default = snapshot.expression_summary_default_cube is not None and compare is None
        for dim in criteria.dict():
            if len(criteria.dict()[dim]) > 0 and depluralize(dim) in expression_summary_non_indexed_dims:
                default = False
                break

        expression_summary = q.expression_summary_default(criteria) if default else q.expression_summary(criteria)

        cell_counts = q.cell_counts(criteria)

        group_by_terms = ["tissue_ontology_term_id", "cell_type_ontology_term_id", compare] if compare else None

        dot_plot_matrix_df, cell_counts_cell_type_agg = get_dot_plot_data(
            expression_summary, cell_counts, group_by_terms
        )
        if is_rollup:
            dot_plot_matrix_df, cell_counts_cell_type_agg = rollup(dot_plot_matrix_df, cell_counts_cell_type_agg)

        response = jsonify(
            dict(
                snapshot_id=snapshot.snapshot_identifier,
                expression_summary=build_expression_summary(dot_plot_matrix_df, compare),
                term_id_labels=dict(
                    genes=build_gene_id_label_mapping(criteria.gene_ontology_term_ids),
                    cell_types=build_ordered_cell_types_by_tissue(
                        cell_counts,
                        cell_counts_cell_type_agg.T,
                        snapshot.cell_type_orderings,
                        compare,
                        group_by_terms,
                    ),
                ),
            )
        )
    return response


def filters():
    request = connexion.request.json
    criteria = WmgFiltersQueryCriteria(**request["filter"])

    with ServerTiming.time("calculate filters and build response"):
        snapshot: WmgSnapshot = load_snapshot()
        response_filter_dims_values = build_filter_dims_values(criteria, snapshot)
        response = jsonify(
            dict(
                snapshot_id=snapshot.snapshot_identifier,
                filter_dims=response_filter_dims_values,
            )
        )
    return response


def markers():
    request = connexion.request.json
    cell_type = request["celltype"]
    tissue = request["tissue"]
    organism = request["organism"]
    n_markers = request["n_markers"]
    test = request["test"]
    snapshot: WmgSnapshot = load_snapshot()

    criteria = MarkerGeneQueryCriteria(
        tissue_ontology_term_id=tissue,
        organism_ontology_term_id=organism,
        cell_type_ontology_term_id=cell_type,
    )
    q = WmgQuery(snapshot)
    df = q.marker_genes(criteria)
    marker_genes = retrieve_top_n_markers(df, test, n_markers)
    return jsonify(
        dict(
            snapshot_id=snapshot.snapshot_identifier,
            marker_genes=marker_genes,
        )
    )


def fetch_datasets_metadata(snapshot: WmgSnapshot, dataset_ids: Iterable[str]) -> List[Dict]:
    return [
        snapshot.dataset_dict.get(dataset_id, dict(id=dataset_id, label="", collection_id="", collection_label=""))
        for dataset_id in dataset_ids
    ]


def find_dimension_id_from_compare(compare: str) -> str:
    if compare == "sex":
        return "sex_ontology_term_id"
    elif compare == "self_reported_ethnicity":
        return "self_reported_ethnicity_ontology_term_id"
    elif compare == "disease":
        return "disease_ontology_term_id"
    else:
        return None


def is_criteria_empty(criteria: WmgFiltersQueryCriteria) -> bool:
    criteria = criteria.dict()
    for key in criteria:
        if key != "organism_ontology_term_id":
            if isinstance(criteria[key], list):
                if len(criteria[key]) > 0:
                    return False
            else:
                if criteria[key] != "":
                    return False
    return True


def build_filter_dims_values(criteria: WmgFiltersQueryCriteria, snapshot: WmgSnapshot) -> Dict:
    dims = {
        "dataset_id": "",
        "disease_ontology_term_id": "",
        "sex_ontology_term_id": "",
        "development_stage_ontology_term_id": "",
        "self_reported_ethnicity_ontology_term_id": "",
        "tissue_ontology_term_id": "",
    }
    for dim in dims:
        dims[dim] = (
            find_all_dim_option_values(snapshot, dim)
            if is_criteria_empty(criteria)
            else find_dim_option_values(criteria, snapshot, dim)
        )

    response_filter_dims_values = dict(
        datasets=fetch_datasets_metadata(snapshot, dims["dataset_id"]),
        disease_terms=build_ontology_term_id_label_mapping(dims["disease_ontology_term_id"]),
        sex_terms=build_ontology_term_id_label_mapping(dims["sex_ontology_term_id"]),
        development_stage_terms=build_ontology_term_id_label_mapping(dims["development_stage_ontology_term_id"]),
        self_reported_ethnicity_terms=build_ontology_term_id_label_mapping(
            dims["self_reported_ethnicity_ontology_term_id"]
        ),
        tissue_terms=build_ontology_term_id_label_mapping(dims["tissue_ontology_term_id"]),
    )

    return response_filter_dims_values


def build_expression_summary(query_result: DataFrame, compare: str) -> dict:
    # Create nested dicts with gene_ontology_term_id, tissue_ontology_term_id keys, cell_type_ontology_term_id respectively
    structured_result: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(dict))
    )

    # Populate aggregated gene expressions
    query_result_agg = query_result.groupby(
        ["gene_ontology_term_id", "tissue_ontology_term_id", "cell_type_ontology_term_id"], as_index=False
    ).agg({"nnz": "sum", "sum": "sum", "n_cells_cell_type": "sum", "n_cells_tissue": "first"})

    for i in range(query_result_agg.shape[0]):
        row = query_result_agg.iloc[i]
        structured_result[row.gene_ontology_term_id][row.tissue_ontology_term_id][row.cell_type_ontology_term_id][
            "aggregated"
        ] = dict(
            n=int(row["nnz"]),
            me=float(row["sum"] / row["nnz"]),
            pc=float(row["nnz"] / row["n_cells_cell_type"]),
            tpc=float(row["nnz"] / row["n_cells_tissue"]),
        )

    # Populate compare filter gene expressions
    if compare:
        for i in range(query_result.shape[0]):
            row = query_result.iloc[i]
            structured_result[row.gene_ontology_term_id][row.tissue_ontology_term_id][row.cell_type_ontology_term_id][
                row[compare]
            ] = dict(
                n=int(row["nnz"]),
                me=float(row["sum"] / row["nnz"]),
                pc=float(row["nnz"] / row["n_cells_cell_type"]),
                tpc=float(row["nnz"] / row["n_cells_tissue"]),
            )

    return structured_result


def agg_cell_type_counts(cell_counts: DataFrame, group_by_terms: List[str] = None) -> DataFrame:
    # Aggregate cube data by tissue, cell type
    if group_by_terms is None:
        group_by_terms = DEFAULT_GROUP_BY_TERMS
    cell_counts_cell_type_agg = cell_counts.groupby(group_by_terms, as_index=True).sum(numeric_only=True)
    cell_counts_cell_type_agg.rename(columns={"n_total_cells": "n_cells_cell_type"}, inplace=True)
    return cell_counts_cell_type_agg


def agg_tissue_counts(cell_counts: DataFrame) -> DataFrame:
    # Aggregate cube data by tissue
    cell_counts_tissue_agg = cell_counts.groupby(["tissue_ontology_term_id"], as_index=True).sum(numeric_only=True)
    cell_counts_tissue_agg.rename(columns={"n_total_cells": "n_cells_tissue"}, inplace=True)
    return cell_counts_tissue_agg


def get_dot_plot_data(
    query_result: DataFrame,
    cell_counts: DataFrame,
    group_by_terms: List[str] = None,
) -> Tuple[DataFrame, DataFrame]:
    if group_by_terms is None:
        group_by_terms = DEFAULT_GROUP_BY_TERMS
    # Get the dot plot matrix dataframe and aggregated cell counts per cell type
    cell_counts_cell_type_agg = agg_cell_type_counts(cell_counts, group_by_terms)
    cell_counts_tissue_agg = agg_tissue_counts(cell_counts)
    dot_plot_matrix_df = build_dot_plot_matrix(
        query_result, cell_counts_cell_type_agg, cell_counts_tissue_agg, group_by_terms
    )
    return dot_plot_matrix_df, cell_counts_cell_type_agg


def rollup(dot_plot_matrix_df, cell_counts_cell_type_agg) -> Tuple[DataFrame, DataFrame]:
    # Roll up numeric columns in the input dataframes
    ignore_cols = ["n_cells_tissue"]

    if dot_plot_matrix_df.shape[0] > 0:
        dot_plot_matrix_df = rollup_across_cell_type_descendants(dot_plot_matrix_df, ignore_cols=ignore_cols)

    if cell_counts_cell_type_agg.shape[0] > 0:
        # make the cell counts dataframe tidy
        for col in cell_counts_cell_type_agg.index.names:
            cell_counts_cell_type_agg[col] = cell_counts_cell_type_agg.index.get_level_values(col)
        cell_counts_cell_type_agg = rollup_across_cell_type_descendants(
            cell_counts_cell_type_agg, ignore_cols=ignore_cols
        )

        # clean up columns that were added to the dataframe to make it tidy
        cell_counts_cell_type_agg.drop(columns=cell_counts_cell_type_agg.index.names, inplace=True)
    return dot_plot_matrix_df, cell_counts_cell_type_agg


def build_dot_plot_matrix(
    query_result: DataFrame,
    cell_counts_cell_type_agg: DataFrame,
    cell_counts_tissue_agg: DataFrame,
    group_by_terms: List[str] = None,
) -> DataFrame:
    if group_by_terms is None:
        group_by_terms = DEFAULT_GROUP_BY_TERMS

    # Aggregate cube data by gene, tissue, cell type
    expr_summary_agg = query_result.groupby(["gene_ontology_term_id"] + group_by_terms, as_index=False).sum(
        numeric_only=True
    )
    return expr_summary_agg.join(cell_counts_cell_type_agg, on=group_by_terms, how="left").join(
        cell_counts_tissue_agg, on=["tissue_ontology_term_id"], how="left"
    )


def build_gene_id_label_mapping(gene_ontology_term_ids: List[str]) -> List[dict]:
    return [
        {gene_ontology_term_id: gene_term_label(gene_ontology_term_id)}
        for gene_ontology_term_id in gene_ontology_term_ids
    ]


def build_ontology_term_id_label_mapping(ontology_term_ids: Iterable[str]) -> List[dict]:
    return [{ontology_term_id: ontology_term_label(ontology_term_id)} for ontology_term_id in ontology_term_ids]


# getting only cell type metadata, no genes
def build_ordered_cell_types_by_tissue(
    cell_counts: DataFrame,
    cell_counts_cell_type_agg_T: DataFrame,
    cell_type_orderings: DataFrame,
    compare: str,
    group_by_terms: List[str] = None,
) -> Dict[str, Dict[str, Dict[str, Any]]]:
    if group_by_terms is None:
        group_by_terms = DEFAULT_GROUP_BY_TERMS

    distinct_tissues_cell_types: DataFrame = cell_counts.groupby(group_by_terms, as_index=False).first()[
        group_by_terms + ["n_total_cells"]
    ]

    # building order for cell types for FE to use
    cell_type_orderings["order"] = range(cell_type_orderings.shape[0])

    # make a multi index
    cell_type_orderings = cell_type_orderings.groupby(["tissue_ontology_term_id", "cell_type_ontology_term_id"]).first()

    indexer = list(
        zip(
            distinct_tissues_cell_types["tissue_ontology_term_id"],
            distinct_tissues_cell_types["cell_type_ontology_term_id"],
        )
    )
    indexer_bool_filter = []
    indexer_filter = []
    for index in indexer:
        indexer_bool_filter.append(index in cell_type_orderings.index)
        if index in cell_type_orderings.index:
            indexer_filter.append(index)

    joined = distinct_tissues_cell_types[indexer_bool_filter]

    for column in cell_type_orderings:
        joined[column] = list(cell_type_orderings[column][indexer])

    # Remove cell types without counts
    joined = joined[joined["n_total_cells"].notnull()]

    # Create nested dicts with tissue_ontology_term_id keys, cell_type_ontology_term_id respectively
    structured_result: Dict[str, Dict[str, Dict[str, Any]]] = defaultdict(lambda: defaultdict(dict))

    # Populate aggregated gene expressions
    joined_agg = joined.groupby(["tissue_ontology_term_id", "cell_type_ontology_term_id"], as_index=False).agg(
        {"n_total_cells": "sum", "depth": "first", "order": "first"}
    )

    agg = cell_counts_cell_type_agg_T.T.groupby(["tissue_ontology_term_id", "cell_type_ontology_term_id"]).sum().T

    for i in range(joined_agg.shape[0]):
        row = joined_agg.iloc[i]
        structured_result[row.tissue_ontology_term_id][row.cell_type_ontology_term_id]["aggregated"] = {
            "cell_type_ontology_term_id": row.cell_type_ontology_term_id,
            "name": ontology_term_label(row.cell_type_ontology_term_id),
            "total_count": int(agg[row.tissue_ontology_term_id][row.cell_type_ontology_term_id]["n_cells_cell_type"]),
            "order": int(row.order),
        }

    # Populate compare filter gene expressions
    if compare:
        for i in range(joined.shape[0]):
            row = joined.iloc[i]
            id_to_label = build_ontology_term_id_label_mapping([row[compare]])[0]
            structured_result[row.tissue_ontology_term_id][row.cell_type_ontology_term_id][row[compare]] = {
                "cell_type_ontology_term_id": row.cell_type_ontology_term_id,
                "name": id_to_label.pop(row[compare]),
                "total_count": int(
                    cell_counts_cell_type_agg_T[row.tissue_ontology_term_id][row.cell_type_ontology_term_id][
                        row[compare]
                    ]["n_cells_cell_type"]
                ),
                "order": int(row.order),
            }

    return structured_result
