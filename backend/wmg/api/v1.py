from collections import defaultdict
from typing import Dict, List, Any, Iterable, Tuple
from math import isnan
from backend.layers.business.business import BusinessLogic
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.common.entities import DatasetId

import connexion
from flask import jsonify
from pandas import DataFrame

from backend.wmg.data.ontology_labels import ontology_term_label, gene_term_label
from backend.wmg.data.query import (
    WmgQuery,
    WmgQueryCriteria,
    MarkerGeneQueryCriteria,
    retrieve_top_n_markers,
)
from backend.wmg.data.snapshot import load_snapshot, WmgSnapshot
from backend.wmg.api.rollup import rollup_across_cell_type_descendants

# TODO: add cache directives: no-cache (i.e. revalidate); impl etag
#  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data
#  -portal/2132


def primary_filter_dimensions():
    snapshot: WmgSnapshot = load_snapshot()
    return jsonify(snapshot.primary_filter_dimensions)


def query():
    request = connexion.request.json
    is_rollup = request.get("is_rollup", True)
    include_filter_dims = request.get("include_filter_dims", False)

    criteria = WmgQueryCriteria(**request["filter"])

    snapshot: WmgSnapshot = load_snapshot()
    q = WmgQuery(snapshot)

    expression_summary = q.expression_summary(criteria)
    cell_counts = q.cell_counts(criteria)
    dot_plot_matrix_df, cell_counts_cell_type_agg = get_dot_plot_data(expression_summary, cell_counts, is_rollup)

    response_filter_dims_values = (
        build_filter_dims_values(criteria, snapshot, expression_summary) if include_filter_dims else {}
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


_business_logic = None


def get_business_logic():
    """
    Returns an instance of the business logic handler. Use this to interrogate the database
    """
    global _business_logic
    if not _business_logic:
        _business_logic = BusinessLogic(DatabaseProvider(), None, None, None, None)
    return _business_logic


# TODO: Read this from generated data artifact instead of DB.
#  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data
#  -portal/2086. This code is without a unit test, but we are intending to replace it.
def fetch_datasets_metadata(dataset_ids: Iterable[str]) -> List[Dict]:
    # We return a DTO because the db entity can't access its attributes after the db session ends,
    # and we want to keep session management out of the calling method

    business_logic = get_business_logic()

    def get_dataset(dataset_id_):

        dataset = business_logic.get_dataset_version_from_canonical(DatasetId(dataset_id_))
        if dataset is None:
            # Handle possible missing dataset due to db state evolving past wmg snapshot
            return dict(id=dataset_id_, label="", collection_id="", collection_label="")
        collection = business_logic.get_collection_version_from_canonical(dataset.collection_id)
        if collection is None:
            # Handle possible missing dataset due to db state evolving past wmg snapshot
            return dict(id=dataset_id_, label="", collection_id="", collection_label="")

        return dict(
            id=dataset.dataset_id.id,
            label=dataset.metadata.name,
            collection_id=collection.collection_id.id,
            collection_label=collection.metadata.name,
        )

    return [get_dataset(dataset_id) for dataset_id in dataset_ids]


def find_dim_option_values(criteria: Dict, snapshot: WmgSnapshot, dimension: str) -> set:
    """Find values for the specified dimension that satisfy the given filtering criteria,
    ignoring any criteria specified for the given dimension."""
    filter_options_criteria = criteria.copy(update={dimension + "s": []}, deep=True)
    # todo can we query cell_counts for a performance gain?
    q = WmgQuery(snapshot)
    query_result = q.expression_summary(filter_options_criteria)
    filter_dims = query_result.groupby(dimension).groups.keys()
    return filter_dims


def build_filter_dims_values(criteria: WmgQueryCriteria, snapshot: WmgSnapshot, expression_summary: DataFrame) -> Dict:
    dims = {
        "dataset_id": "",
        "disease_ontology_term_id": "",
        "sex_ontology_term_id": "",
        "development_stage_ontology_term_id": "",
        "self_reported_ethnicity_ontology_term_id": "",
    }
    for dim in dims:
        if len(criteria.dict()[dim + "s"]) == 0:
            dims[dim] = expression_summary.groupby(dim).groups.keys()
        else:
            dims[dim] = find_dim_option_values(criteria, snapshot, dim)

    response_filter_dims_values = dict(
        datasets=fetch_datasets_metadata(dims["dataset_id"]),
        disease_terms=build_ontology_term_id_label_mapping(dims["disease_ontology_term_id"]),
        sex_terms=build_ontology_term_id_label_mapping(dims["sex_ontology_term_id"]),
        development_stage_terms=build_ontology_term_id_label_mapping(dims["development_stage_ontology_term_id"]),
        self_reported_ethnicity_terms=build_ontology_term_id_label_mapping(
            dims["self_reported_ethnicity_ontology_term_id"]
        ),
    )

    return response_filter_dims_values


def build_expression_summary(query_result: DataFrame) -> dict:
    # Create nested dicts with gene_ontology_term_id, tissue_ontology_term_id keys, respectively
    # is_rollup is a flag to indicate whether the expressions should be rolled up or not
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


def get_dot_plot_data(query_result: DataFrame, cell_counts: DataFrame, is_rollup: bool) -> Tuple[DataFrame, DataFrame]:
    # Get the dot plot matrix dataframe and aggregated cell counts per cell type
    cell_counts_cell_type_agg = agg_cell_type_counts(cell_counts)
    cell_counts_tissue_agg = agg_tissue_counts(cell_counts)
    dot_plot_matrix_df = build_dot_plot_matrix(query_result, cell_counts_cell_type_agg, cell_counts_tissue_agg)

    if is_rollup:
        dot_plot_matrix_df = rollup_across_cell_type_descendants(dot_plot_matrix_df)
        # make the cell counts dataframe tidy
        for col in cell_counts_cell_type_agg.index.names:
            cell_counts_cell_type_agg[col] = cell_counts_cell_type_agg.index.get_level_values(col)
        cell_counts_cell_type_agg = rollup_across_cell_type_descendants(cell_counts_cell_type_agg)
        # clean up columns that were added to the dataframe to make it tidy
        for col in cell_counts_cell_type_agg.index.names:
            del cell_counts_cell_type_agg[col]

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
