from collections import defaultdict
from typing import Dict, List, Any, Iterable
from math import isnan

import connexion
from flask import jsonify
from pandas import DataFrame

from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.db_session import db_session_manager
from backend.wmg.data.ontology_labels import ontology_term_label, gene_term_label
from backend.wmg.data.query import WmgQuery, WmgQueryCriteria
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
    dot_plot_matrix_df = build_dot_plot_matrix(expression_summary, cell_counts)

    include_filter_dims = request.get("include_filter_dims", False)
    response_filter_dims_values = (
        build_filter_dims_values(criteria, snapshot.filter_relationships) if include_filter_dims else {}
    )
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


def find_dim_option_values(criteria: Dict, htable: Dict, dimension: str) -> set:
    """Find values for the specified dimension that satisfy the given filtering criteria,
    ignoring any criteria specified for the given dimension."""

    dimension = dimension[:-1] if dimension[-1] == "s" else dimension
    supersets = []
    for key in criteria:
        attrs = criteria[key]
        key = key[:-1] if key[-1] == "s" else key
        if key != dimension:
            if isinstance(attrs, list):
                if len(attrs) > 0:
                    vals = [key + "__" + val for val in attrs]
                    sets = [set([x for x in htable[v] if dimension in x]) for v in vals]
                    supersets.append(sets[0].union(*sets[1:]))
            else:
                if attrs != "":
                    supersets.append(set([x for x in htable[key + "__" + attrs] if dimension in x]))

    if len(supersets) > 1:
        valid_options = supersets[0].intersection(*supersets[1:])
    else:
        valid_options = supersets[0]

    final_options = []
    for v in valid_options:
        loop_back_options = htable[v]
        if len(set(loop_back_options).intersection(*supersets)) > 0:
            final_options.append(v)

    return [i.split("__")[1] for i in final_options]


def build_filter_dims_values(criteria: WmgQueryCriteria, htable: Dict) -> Dict:
    dims = {
        "dataset_id": "",
        "disease_ontology_term_id": "",
        "sex_ontology_term_id": "",
        "development_stage_ontology_term_id": "",
        "ethnicity_ontology_term_id": "",
    }

    criteria = dict(criteria)
    del criteria["gene_ontology_term_ids"]

    for dim in dims:
        dims[dim] = find_dim_option_values(criteria, htable, dim)

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
                "total_count": row.n_total_cells,
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
