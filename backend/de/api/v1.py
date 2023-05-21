from typing import Dict, Iterable, List

import connexion
import numpy as np
import pandas as pd
from flask import jsonify
from server_timing import Timing as ServerTiming

from backend.wmg.data.calculate_markers import _run_ttest
from backend.wmg.data.ontology_labels import ontology_term_label
from backend.wmg.data.query import (
    FmgQueryCriteria,
    WmgFiltersQueryCriteria,
    WmgQuery
)
from backend.wmg.data.schemas.cube_schema import expression_summary_non_indexed_dims
from backend.wmg.data.snapshot import WmgSnapshot, load_snapshot
from backend.wmg.data.utils import depluralize, find_all_dim_option_values, find_dim_option_values, to_dict


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

def differentialExpression():
    request = connexion.request.json

    queryGroup1Filters = request["queryGroup1Filters"]
    queryGroup2Filters = request["queryGroup2Filters"]

    criteria1 = FmgQueryCriteria(**queryGroup1Filters)
    criteria2 = FmgQueryCriteria(**queryGroup2Filters)

    snapshot: WmgSnapshot = load_snapshot()

    q = WmgQuery(snapshot)

    with ServerTiming.time("run differential expression"):
        results1, results2, success_code = run_differential_expression(q, criteria1, criteria2)

    return jsonify(
        dict(
            snapshot_id=snapshot.snapshot_identifier,
            differentialExpressionResults1=results1,
            differentialExpressionResults2=results2,
            successCode=success_code,
        )
    )


def fetch_datasets_metadata(snapshot: WmgSnapshot, dataset_ids: Iterable[str]) -> List[Dict]:
    return [
        snapshot.dataset_dict.get(dataset_id, dict(id=dataset_id, label="", collection_id="", collection_label=""))
        for dataset_id in dataset_ids
    ]

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
        cell_type_terms=build_ontology_term_id_label_mapping(dims["cell_type_ontology_term_id"]),
    )

    return response_filter_dims_values


def build_ontology_term_id_label_mapping(ontology_term_ids: Iterable[str]) -> List[dict]:
    return [{ontology_term_id: ontology_term_label(ontology_term_id)} for ontology_term_id in ontology_term_ids]

def should_use_default_cube(criteria):
    default = True
    for dim in criteria.dict():
        if len(criteria.dict()[dim]) > 0 and depluralize(dim) in expression_summary_non_indexed_dims:
            default = False
            break
    return default


def de_get_expression_summary(q, criteria):
    if should_use_default_cube(criteria):
        es_agg = q.expression_summary_default(criteria).groupby("gene_ontology_term_id").sum(numeric_only=True)
    else:
        es_agg = q.expression_summary(criteria).groupby("gene_ontology_term_id").sum(numeric_only=True)

    return es_agg


def run_differential_expression(q, criteria1, criteria2, pval_thr=1e-5, threshold=2000):
    cell_counts1 = q.cell_counts(criteria1)
    cell_counts2 = q.cell_counts(criteria2)

    n_cells1 = cell_counts1["n_total_cells"].sum()
    n_cells2 = cell_counts2["n_total_cells"].sum()

    skip1 = not should_use_default_cube(criteria1) and cell_counts1.shape[0] > threshold
    skip2 = not should_use_default_cube(criteria2) and cell_counts2.shape[0] > threshold

    skip = skip1 or skip2

    success_code = 0
    if skip1 and skip2:
        success_code = 3
    elif skip1:
        success_code = 1
    elif skip2:
        success_code = 2

    if skip:
        return [], [], success_code

    es_agg1 = de_get_expression_summary(q, criteria1)
    es_agg2 = de_get_expression_summary(q, criteria2)

    genes = list(set(list(es_agg1.index) + list(es_agg2.index)))

    genes_indexer = pd.Series(index=genes, data=np.arange(len(genes)))

    sums1 = np.zeros(len(genes))
    sqsums1 = np.zeros(len(genes))

    sums2 = np.zeros(len(genes))
    sqsums2 = np.zeros(len(genes))

    sums1[genes_indexer[es_agg1.index]] = es_agg1["sum"].values
    sqsums1[genes_indexer[es_agg1.index]] = es_agg1["sqsum"].values

    sums2[genes_indexer[es_agg2.index]] = es_agg2["sum"].values
    sqsums2[genes_indexer[es_agg2.index]] = es_agg2["sqsum"].values

    pvals, effects = _run_ttest(sums1, sqsums1, n_cells1, sums2, sqsums2, n_cells2)
    de_genes = np.array(genes)[np.argsort(-effects)]
    p = pvals[np.argsort(-effects)]
    effects = effects[np.argsort(-effects)]

    statistics1 = []
    for i in range(len(p)):
        pi = p[i]
        ei = abs(effects[i])
        if ei is not np.nan and pi is not np.nan and pi < pval_thr:
            statistics1.append({"gene_ontology_term_id": de_genes[i], "p_value": pi, "effect_size": ei})
            if len(statistics1) >= 250:
                break

    pvals, effects = _run_ttest(sums2, sqsums2, n_cells2, sums1, sqsums1, n_cells1)
    de_genes = np.array(genes)[np.argsort(-effects)]
    p = pvals[np.argsort(-effects)]
    effects = effects[np.argsort(-effects)]

    statistics2 = []
    for i in range(len(p)):
        pi = p[i]
        ei = abs(effects[i])
        if ei is not np.nan and pi is not np.nan and pi < pval_thr:
            statistics2.append({"gene_ontology_term_id": de_genes[i], "p_value": pi, "effect_size": ei})
            if len(statistics2) >= 250:
                break
    return statistics1, statistics2, success_code
