from collections import defaultdict
from typing import Dict, List

import connexion
from flask import jsonify
from pandas import DataFrame

from backend.corpora.common.corpora_orm import DbDataset
from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.db_session import db_session_manager
from backend.wmg.data.cube import load_cube
from backend.wmg.data.query import (
    build_gene_id_label_mapping,
    build_ontology_term_id_label_mapping,
    WmgQuery,
    WmgQueryCriteria,
    build_dot_plot_matrix,
)
from backend.wmg.data.schema import cube_non_indexed_dims


# TODO: add cache directives: no-cache (i.e. revalidate); impl etag
def primary_filter_dimensions():
    cube, snapshot_identifier = load_cube()
    qry = WmgQuery(cube)

    # gene terms are grouped by organism, and represented as a nested lists in dict, keyed by organism
    organism_gene_ids: dict[str, List[str]] = qry.list_grouped_primary_filter_dimensions_term_ids(
        "gene_ontology_term_id", group_by_dim="organism_ontology_term_id"
    )
    organism_gene_terms = {
        organism_term_id: build_gene_id_label_mapping(gene_term_ids)
        for organism_term_id, gene_term_ids in organism_gene_ids.items()
    }

    result = dict(
        snapshot_id=snapshot_identifier,
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

    cube, snapshot_identifier = load_cube()
    query_result = WmgQuery(cube).expression_summary(criteria)
    dot_plot_matrix_df = build_dot_plot_matrix(query_result)
    all_filter_dims_values = extract_filter_dims_values(query_result)

    include_filter_dims = request.get("include_filter_dims", False)
    response_filter_dims_values = build_filter_dims_values(all_filter_dims_values) if include_filter_dims else {}

    cell_type_term_ids = all_filter_dims_values["cell_type_ontology_term_id"]

    return jsonify(
        dict(
            snapshot_id=snapshot_identifier,
            expression_summary=build_expression_summary(dot_plot_matrix_df),
            term_id_labels=dict(
                genes=build_gene_id_label_mapping(criteria.gene_ontology_term_ids),
                cell_types=build_ontology_term_id_label_mapping(cell_type_term_ids),
            ),
            filter_dims=response_filter_dims_values,
        )
    )


# TODO: It would be reasonable to fetch datasets via a REST API call to Data Portal, rather than coupling to the
#  portal's db layer, considering that this is the only occasion where the backend.wmg.api package introduces a
#  dependency on the db. There is no appropriate API call at this time and the lower of making an API call would
#  have to be considered as well.
def fetch_datasets(dataset_ids: List[str]) -> List[Dict]:
    # We return a DTO because the db entity can't access its attributes after the db session ends,
    # and we want to keep session management out of the calling method
    def result(dataset: DbDataset):
        return dict(
            id=dataset.id, name=dataset.name, collection=dict(id=dataset.collection.id, name=dataset.collection.name)
        )

    with db_session_manager() as session:
        return [result(Dataset.get(session, dataset_id).db_object) for dataset_id in dataset_ids]


def build_datasets(dataset_ids: List[str]) -> List[Dict]:
    datasets = fetch_datasets(dataset_ids)
    return [
        dict(
            id=dataset["id"],
            label=dataset["name"],
            collection_id=dataset["collection"]["id"],
            collection_label=dataset["collection"]["name"],
        )
        for dataset in datasets
    ]


def build_filter_dims_values(all_filter_dims_values):
    response_filter_dims_values = dict(
        datasets=build_datasets(all_filter_dims_values["dataset_id"]),
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
    structured_result = defaultdict(lambda: defaultdict(list))
    for group_by_key, cell_type_stats in query_result.to_dict("index").items():
        gene_ontology_term_id, tissue_ontology_term_id, cell_type_ontology_term_id = [s for s in group_by_key]
        structured_result[gene_ontology_term_id][tissue_ontology_term_id].append(
            dict(
                id=cell_type_ontology_term_id,
                n=cell_type_stats["nnz"],
                me=cell_type_stats["sum"] / cell_type_stats["n_cells"],
                # TODO
                pc=0.0,
                # TODO
                tpc=0.0,
            )
        )
    return structured_result


def extract_filter_dims_values(query_result: DataFrame) -> dict:
    """
    Return unique values for each dimension in the specified query result
    """
    return {
        col: query_result.groupby(col).groups.keys() for col in set(query_result.columns) & set(cube_non_indexed_dims)
    }
