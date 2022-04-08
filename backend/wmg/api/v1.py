from collections import defaultdict
from functools import cache
from typing import Dict, List, Any, Iterable

import connexion
from flask import jsonify
from pandas import DataFrame

from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.db_session import db_session_manager
from backend.wmg.data.ontology_labels import ontology_term_label, gene_term_label
from backend.wmg.data.query import (
    WmgQuery,
    WmgQueryCriteria,
)
from backend.wmg.data.schemas.cube_schema import cube_non_indexed_dims
from backend.wmg.data.snapshot import load_snapshot, WmgSnapshot


# TODO: add cache directives: no-cache (i.e. revalidate); impl etag
#  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data
#  -portal/2132


def primary_filter_dimensions():
    snapshot: WmgSnapshot = load_snapshot()
    qry = WmgQuery(snapshot)

    # gene terms are grouped by organism, and represented as a nested lists in dict, keyed by organism
    organism_gene_ids: dict[str, List[str]] = qry.list_grouped_primary_filter_dimensions_term_ids(
        "gene_ontology_term_id", group_by_dim="organism_ontology_term_id"
    )
    organism_gene_terms = {
        organism_term_id: build_gene_id_label_mapping(gene_term_ids)
        for organism_term_id, gene_term_ids in organism_gene_ids.items()
    }

    result = dict(
        snapshot_id=snapshot.snapshot_identifier,
        organism_terms=build_ontology_term_id_label_mapping(
            qry.list_primary_filter_dimension_term_ids("organism_ontology_term_id")
        ),
        tissue_terms=build_ontology_term_id_label_mapping(
            qry.list_primary_filter_dimension_term_ids("tissue_ontology_term_id")
        ),
        gene_terms=organism_gene_terms,
    )
    return jsonify(result)


def query():
    request = connexion.request.json
    criteria = WmgQueryCriteria(**request["filter"])

    snapshot: WmgSnapshot = load_snapshot()
    query = WmgQuery(snapshot)
    expression_summary = query.expression_summary(criteria)
    cell_counts = query.cell_counts(criteria)
    dot_plot_matrix_df = build_dot_plot_matrix(expression_summary, cell_counts)
    all_filter_dims_values = extract_filter_dims_values(expression_summary)

    include_filter_dims = request.get("include_filter_dims", False)
    response_filter_dims_values = build_filter_dims_values(all_filter_dims_values) if include_filter_dims else {}

    return jsonify(
        dict(
            snapshot_id=snapshot.snapshot_identifier,
            expression_summary=build_expression_summary(dot_plot_matrix_df),
            term_id_labels=dict(
                genes=build_gene_id_label_mapping(criteria.gene_ontology_term_ids),
                cell_types=build_ordered_cell_types_by_tissue(cell_counts, snapshot.cell_type_orderings),
            ),
            filter_dims=response_filter_dims_values,
        )
    )


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


def build_filter_dims_values(all_filter_dims_values: Dict[str, Iterable[str]]):
    response_filter_dims_values = dict(
        datasets=fetch_datasets_metadata(all_filter_dims_values["dataset_id"]),
        disease_terms=build_ontology_term_id_label_mapping(all_filter_dims_values["disease_ontology_term_id"]),
        sex_terms=build_ontology_term_id_label_mapping(all_filter_dims_values["sex_ontology_term_id"]),
        development_stage_terms=build_ontology_term_id_label_mapping(
            all_filter_dims_values["development_stage_ontology_term_id"]
        ),
        ethnicity_terms=build_ontology_term_id_label_mapping(all_filter_dims_values["ethnicity_ontology_term_id"]),
        # excluded per product requirements, but keeping in, commented-out, to reduce future head-scratching
        # assay_ontology_terms=build_ontology_term_id_label_mapping(all_filter_dims_values["assay_ontology_term_id"]),
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


def extract_filter_dims_values(query_result: DataFrame) -> Dict[str, set]:
    """
    Return unique values for each dimension in the specified query result
    """
    dims: set = set(query_result.columns) & set(cube_non_indexed_dims)
    dim_uniq_values: Dict[str, set] = {dim: query_result.groupby(dim).groups.keys() for dim in dims}
    return dim_uniq_values


def build_dot_plot_matrix(query_result: DataFrame, cell_counts: DataFrame) -> DataFrame:
    # Aggregate cube data by gene, tissue, cell type
    expr_summary_agg = query_result.groupby(
        ["gene_ontology_term_id", "tissue_ontology_term_id", "cell_type_ontology_term_id"], as_index=False
    ).sum()

    cell_counts_cell_type_agg = cell_counts.groupby(
        ["tissue_ontology_term_id", "cell_type_ontology_term_id"], as_index=True
    ).sum()
    cell_counts_cell_type_agg.rename(columns={"n_total_cells": "n_cells_cell_type"}, inplace=True)

    cell_counts_tissue_agg = cell_counts.groupby(["tissue_ontology_term_id"], as_index=True).sum()
    cell_counts_tissue_agg.rename(columns={"n_total_cells": "n_cells_tissue"}, inplace=True)

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
    cell_counts: DataFrame, cell_type_orderings: DataFrame
) -> Dict[str, List[Dict[str, str]]]:
    distinct_tissues_cell_types: DataFrame = cell_counts.groupby(
        ["tissue_ontology_term_id", "cell_type_ontology_term_id"], as_index=False
    ).first()[["tissue_ontology_term_id", "cell_type_ontology_term_id"]]

    joined = distinct_tissues_cell_types.merge(
        cell_type_orderings, on=["tissue_ontology_term_id", "cell_type_ontology_term_id"]
    )
    sorted = joined.sort_values(by=["tissue_ontology_term_id", "order"])

    structured_result: Dict[str, List[Dict[str, str]]] = defaultdict(list)

    for row in sorted.itertuples(index=False):
        structured_result[row.tissue_ontology_term_id].append(
            {row.cell_type_ontology_term_id: ontology_term_label(row.cell_type_ontology_term_id)}
        )

    return structured_result
