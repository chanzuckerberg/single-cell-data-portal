from collections import defaultdict
from typing import Any, Dict, Iterable, List

import connexion
from flask import jsonify
from pandas import DataFrame
from server_timing import Timing as ServerTiming

from backend.wmg.api.common.expression_dotplot import get_dot_plot_data
from backend.wmg.api.common.rollup import rollup
from backend.wmg.api.wmg_api_config import WMG_API_FORCE_LOAD_SNAPSHOT_ID, WMG_API_SNAPSHOT_SCHEMA_VERSION
from backend.wmg.data.ontology_labels import gene_term_label, ontology_term_label
from backend.wmg.data.query import (
    MarkerGeneQueryCriteria,
    WmgFiltersQueryCriteria,
    WmgQuery,
    WmgQueryCriteriaV2,
    retrieve_top_n_markers,
)
from backend.wmg.data.schemas.cube_schema import expression_summary_non_indexed_dims
from backend.wmg.data.snapshot import WmgSnapshot, load_snapshot
from backend.wmg.data.utils import depluralize, find_all_dim_option_values, find_dim_option_values

# TODO: add cache directives: no-cache (i.e. revalidate); impl etag
#  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data
#  -portal/2132


def primary_filter_dimensions():
    snapshot: WmgSnapshot = load_snapshot(
        snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION,
        explicit_snapshot_id_to_load=WMG_API_FORCE_LOAD_SNAPSHOT_ID,
    )

    return jsonify(snapshot.primary_filter_dimensions)


def query():
    request = connexion.request.json
    is_rollup = request.get("is_rollup", True)
    compare = request.get("compare", None)

    if compare:
        compare = find_dimension_id_from_compare(compare)

    criteria = WmgQueryCriteriaV2(**request["filter"])

    with ServerTiming.time("query and build response"):
        snapshot: WmgSnapshot = load_snapshot(
            snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION,
            explicit_snapshot_id_to_load=WMG_API_FORCE_LOAD_SNAPSHOT_ID,
        )

        q = WmgQuery(snapshot)
        default = snapshot.expression_summary_default_cube is not None and compare is None
        for dim in criteria.dict():
            if len(criteria.dict()[dim]) > 0 and depluralize(dim) in expression_summary_non_indexed_dims:
                default = False
                break

        expression_summary = (
            q.expression_summary_default(criteria)
            if default
            else q.expression_summary(criteria, compare_dimension=compare)
        )

        cell_counts = q.cell_counts(criteria, compare_dimension=compare)
        if expression_summary.shape[0] > 0 or cell_counts.shape[0] > 0:
            group_by_terms = ["tissue_ontology_term_id", "cell_type_ontology_term_id", compare] if compare else None

            gene_expression_df, cell_counts_grouped_df = get_dot_plot_data(
                expression_summary, cell_counts, group_by_terms
            )
            if is_rollup:
                rolled_gene_expression_df, rolled_cell_counts_grouped_df = rollup(
                    gene_expression_df, cell_counts_grouped_df
                )

            response = jsonify(
                dict(
                    snapshot_id=snapshot.snapshot_identifier,
                    expression_summary=build_expression_summary(gene_expression_df, rolled_gene_expression_df, compare),
                    term_id_labels=dict(
                        genes=build_gene_id_label_mapping(criteria.gene_ontology_term_ids),
                        cell_types=build_ordered_cell_types_by_tissue(
                            rolled_cell_counts_grouped_df, snapshot.cell_type_orderings, compare
                        ),
                    ),
                )
            )
        else:  # no data, return empty json
            response = jsonify(dict(snapshot_id=snapshot.snapshot_identifier, expression_summary={}, term_id_labels={}))
    return response


def filters():
    request = connexion.request.json
    criteria = WmgFiltersQueryCriteria(**request["filter"])

    with ServerTiming.time("calculate filters and build response"):
        snapshot: WmgSnapshot = load_snapshot(
            snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION,
            explicit_snapshot_id_to_load=WMG_API_FORCE_LOAD_SNAPSHOT_ID,
        )

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
    snapshot: WmgSnapshot = load_snapshot(
        snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION,
        explicit_snapshot_id_to_load=WMG_API_FORCE_LOAD_SNAPSHOT_ID,
    )

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
    elif compare == "publication":
        return "publication_citation"
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
        "cell_type_ontology_term_id": "",
        "publication_citation": "",
    }
    for dim in dims:
        dims[dim] = (
            find_all_dim_option_values(snapshot, criteria.organism_ontology_term_id, dim)
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
        publication_citations=dims["publication_citation"],
        cell_type_terms=build_ontology_term_id_label_mapping(dims["cell_type_ontology_term_id"]),
        tissue_terms=build_ontology_term_id_label_mapping(dims["tissue_ontology_term_id"]),
    )

    return response_filter_dims_values


def build_expression_summary(
    unrolled_gene_expression_df: DataFrame, rolled_gene_expression_df: DataFrame, compare: str
) -> dict:
    """
    Compute and build a data structure that contains gene expression summary statistics.

    Parameters
    ----------
    unrolled_gene_expression_df: A dataframe containing unrolled gene expression values for each gene.

    rolled_gene_expression_df: A dataframe containing rolleup gene expression values.

    compare: Optional. The compare dimension to further group the gene expression summary statistics into.

    Returns
    -------
    structured_result : A nested dictionary that contains gene expression summary statistics.
    """
    # Create nested dicts with gene_ontology_term_id, tissue_ontology_term_id keys, cell_type_ontology_term_id respectively
    structured_result: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(dict))
    )

    # Populate gene expressions stats for each (gene, tissue) combination
    #
    # For aggregations that do not contain `cell_type_ontology_term_id` as a component in the compound key
    # used to group rows, we SHOULD NOT USE the gene expression dataframe that contains rolled up values.
    # That is, we should not use `rolled_gene_expression_df` for such aggregations.
    # This is because the roll up is computed over the cell_type ontology inheritance graph,
    # and such aggregations will count the rolled up values many times.
    #
    # For group-by compound keys that do not contain the `cell_type_ontology_term_id`,
    # the unrolled gene expression dataframe containing unaggregated gene expression values should be used.
    # That is, we use `unrolled_gene_expression_df` for aggregation over the group by compound key:
    # (`gene_ontology_term_id`, `tissue_ontology_term_id`)
    gene_expr_grouped_df = unrolled_gene_expression_df.groupby(
        ["gene_ontology_term_id", "tissue_ontology_term_id"], as_index=False
    ).agg({"nnz": "sum", "sum": "sum", "n_cells_tissue": "first"})

    for i in range(gene_expr_grouped_df.shape[0]):
        row = gene_expr_grouped_df.iloc[i]
        structured_result[row.gene_ontology_term_id][row.tissue_ontology_term_id]["tissue_stats"]["aggregated"] = dict(
            n=int(row["nnz"]),
            me=(float(row["sum"] / row["nnz"]) if row["nnz"] else 0.0),
            tpc=(float(row["nnz"] / row["n_cells_tissue"]) if row["n_cells_tissue"] else 0.0),
        )

    # Populate gene expressions stats for each (gene, tissue, cell_type) combination
    #
    # For aggregations that contain `cell_type_ontology_term_id` as a component in the compound key used
    # to group rows, we can use the gene expression dataframe that contains rolled up values.
    # That is, we can use `rolled_gene_expression_df` for such aggregations.
    gene_expr_grouped_df = rolled_gene_expression_df.groupby(
        ["gene_ontology_term_id", "tissue_ontology_term_id", "cell_type_ontology_term_id"], as_index=False
    ).agg({"nnz": "sum", "sum": "sum", "n_cells_cell_type": "sum", "n_cells_tissue": "first"})

    fill_out_structured_dict_aggregated(gene_expr_grouped_df, structured_result)

    # Populate gene expression stats for each (gene, tissue, cell_type, <compare_dimension>) combination
    if compare:
        fill_out_structured_dict_compare(rolled_gene_expression_df, structured_result, compare)

    return structured_result


def fill_out_structured_dict_aggregated(gene_expr_grouped_df, structured_result):
    n = gene_expr_grouped_df["nnz"].astype("int").values
    me = (gene_expr_grouped_df["sum"] / gene_expr_grouped_df["nnz"]).values
    pc = (gene_expr_grouped_df["nnz"] / gene_expr_grouped_df["n_cells_cell_type"]).values
    tpc = (gene_expr_grouped_df["nnz"] / gene_expr_grouped_df["n_cells_tissue"]).values
    genes = gene_expr_grouped_df["gene_ontology_term_id"].values
    tissues = gene_expr_grouped_df["tissue_ontology_term_id"].values
    cell_types = gene_expr_grouped_df["cell_type_ontology_term_id"].values

    for i in range(len(n)):
        structured_result[genes[i]][tissues[i]][cell_types[i]]["aggregated"] = dict(
            n=int(n[i]),
            me=float(me[i]),
            pc=float(pc[i]),
            tpc=float(tpc[i]),
        )


def fill_out_structured_dict_compare(rolled_gene_expression_df, structured_result, compare):
    n = rolled_gene_expression_df["nnz"].astype("int").values
    me = (rolled_gene_expression_df["sum"] / rolled_gene_expression_df["nnz"]).values
    pc = (rolled_gene_expression_df["nnz"] / rolled_gene_expression_df["n_cells_cell_type"]).values
    tpc = (rolled_gene_expression_df["nnz"] / rolled_gene_expression_df["n_cells_tissue"]).values
    genes = rolled_gene_expression_df["gene_ontology_term_id"].values
    tissues = rolled_gene_expression_df["tissue_ontology_term_id"].values
    cell_types = rolled_gene_expression_df["cell_type_ontology_term_id"].values
    compare = rolled_gene_expression_df[compare].values

    for i in range(len(n)):
        structured_result[genes[i]][tissues[i]][cell_types[i]][compare[i]] = dict(
            n=int(n[i]),
            me=float(me[i]),
            pc=float(pc[i]),
            tpc=float(tpc[i]),
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
    cell_counts_cell_type_agg: DataFrame,
    cell_type_orderings: DataFrame,
    compare: str,
) -> Dict[str, Dict[str, Dict[str, Any]]]:

    distinct_tissues_cell_types = cell_counts_cell_type_agg.reset_index().rename(
        columns={"n_cells_cell_type": "n_total_cells"}
    )

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
        joined[column] = list(cell_type_orderings[column][indexer_filter])

    # Remove cell types without counts
    joined = joined[joined["n_total_cells"].notnull()]

    # Create nested dicts with tissue_ontology_term_id keys, cell_type_ontology_term_id respectively
    structured_result: Dict[str, Dict[str, Dict[str, Any]]] = defaultdict(lambda: defaultdict(dict))

    # Populate aggregated cell counts for each tissue
    joined_agg = joined.groupby(["tissue_ontology_term_id"], as_index=False).agg({"n_total_cells": "sum"})

    for i in range(joined_agg.shape[0]):
        row = joined_agg.iloc[i]
        structured_result[row.tissue_ontology_term_id]["tissue_stats"]["aggregated"] = {
            "tissue_ontology_term_id": row.tissue_ontology_term_id,
            "name": ontology_term_label(row.tissue_ontology_term_id),
            "total_count": int(row.n_total_cells),
            "order": -1,
        }

    # Populate aggregated cell counts for each (tissue, cell_type) combination
    joined_agg = joined.groupby(["tissue_ontology_term_id", "cell_type_ontology_term_id"], as_index=False).agg(
        {"order": "first"}
    )

    agg = cell_counts_cell_type_agg.groupby(["tissue_ontology_term_id", "cell_type_ontology_term_id"]).sum().T

    for i in range(joined_agg.shape[0]):
        row = joined_agg.iloc[i]
        structured_result[row.tissue_ontology_term_id][row.cell_type_ontology_term_id]["aggregated"] = {
            "cell_type_ontology_term_id": row.cell_type_ontology_term_id,
            "name": ontology_term_label(row.cell_type_ontology_term_id),
            "total_count": int(agg[row.tissue_ontology_term_id][row.cell_type_ontology_term_id]["n_cells_cell_type"]),
            "order": int(row.order),
        }

    # Populate aggregated cell counts for each (tissue, cell_type, <compare_dimension>) combination
    cell_counts_cell_type_agg_T = cell_counts_cell_type_agg.T
    if compare:
        for i in range(joined.shape[0]):
            row = joined.iloc[i]
            id_to_label = build_ontology_term_id_label_mapping([row[compare]])[0]
            name = id_to_label.pop(row[compare])
            structured_result[row.tissue_ontology_term_id][row.cell_type_ontology_term_id][row[compare]] = {
                "cell_type_ontology_term_id": row.cell_type_ontology_term_id,
                "name": name if name else row[compare],
                "total_count": int(
                    cell_counts_cell_type_agg_T[row.tissue_ontology_term_id][row.cell_type_ontology_term_id][
                        row[compare]
                    ]["n_cells_cell_type"]
                ),
                "order": int(row.order),
            }

    return structured_result
