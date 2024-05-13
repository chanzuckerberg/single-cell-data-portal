import os
from typing import Any, Dict, Iterable, List, Tuple

import connexion
import numpy as np
import pandas as pd
from ddtrace import tracer
from flask import jsonify
from scipy import stats
from server_timing import Timing as ServerTiming

from backend.common.utils.rollup import descendants
from backend.wmg.api.wmg_api_config import (
    WMG_API_FORCE_LOAD_SNAPSHOT_ID,
    WMG_API_READ_FS_CACHED_SNAPSHOT,
    WMG_API_SNAPSHOT_FS_CACHE_ROOT_PATH,
    WMG_API_SNAPSHOT_SCHEMA_VERSION,
)
from backend.wmg.data.ontology_labels import gene_term_label, ontology_term_label
from backend.wmg.data.query import DeQueryCriteria, WmgFiltersQueryCriteria, WmgQuery
from backend.wmg.data.snapshot import WmgSnapshot, load_snapshot
from backend.wmg.data.utils import (
    find_all_dim_option_values,
    find_dim_option_values,
)

DEPLOYMENT_STAGE = os.environ.get("DEPLOYMENT_STAGE", "")
SNAPSHOT_FS_ROOT_PATH = (
    WMG_API_SNAPSHOT_FS_CACHE_ROOT_PATH if (WMG_API_READ_FS_CACHED_SNAPSHOT and DEPLOYMENT_STAGE != "test") else None
)


@tracer.wrap(name="filters", service="wmg-api", resource="filters", span_type="wmg-api")
def filters():
    request = connexion.request.json
    sanitize_api_query_dict(request["filter"])

    criteria = WmgFiltersQueryCriteria(**request["filter"])

    with ServerTiming.time("load snapshot"):
        snapshot: WmgSnapshot = load_snapshot(
            snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION,
            explicit_snapshot_id_to_load=WMG_API_FORCE_LOAD_SNAPSHOT_ID,
            snapshot_fs_root_path=SNAPSHOT_FS_ROOT_PATH,
        )

    with ServerTiming.time("calculate filters and build response"):
        response_filter_dims_values = build_filter_dims_values(criteria, snapshot)
        response = jsonify(
            dict(
                snapshot_id=snapshot.snapshot_identifier,
                filter_dims=response_filter_dims_values,
            )
        )
    return response


def sanitize_api_query_dict(query_dict: Any):
    """
    Remove invalid values in the query dictionary encoding the query API
    request body.

    The assumption is that this function is called at the beginning of the
    API function. This usage also helps mitigate query injection attacks.

    NOTE: This is a destructive operation in that it mutates `query_dict`.

    Parameters
    ----------
    query_dict : json object
        The query dictionary to sanitize.

    Returns
    -------
    None because this function mutates the function argument
    """

    # Sanitize `self_reported_ethnicity_ontology_term_ids` by removing
    # comma-delimited values because WMG does not support filtering and grouping
    # by ethnicity terms that encode mixed ethnicities encoded as a single comma-delimited string
    # value
    if "self_reported_ethnicity_ontology_term_ids" in query_dict:
        ethnicity_term_ids = query_dict["self_reported_ethnicity_ontology_term_ids"]

        ethnicity_term_ids_to_keep = [x for x in ethnicity_term_ids if "," not in x]
        query_dict["self_reported_ethnicity_ontology_term_ids"] = ethnicity_term_ids_to_keep


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


@tracer.wrap(name="build_filter_dims_values", service="wmg-api", resource="filters", span_type="wmg-api")
def build_filter_dims_values(criteria: WmgFiltersQueryCriteria, snapshot: WmgSnapshot) -> Dict:
    dims = {
        "publication_citation": "",
        "disease_ontology_term_id": "",
        "sex_ontology_term_id": "",
        "development_stage_ontology_term_id": "",
        "self_reported_ethnicity_ontology_term_id": "",
        "tissue_ontology_term_id": "",
        "cell_type_ontology_term_id": "",
        "organism_ontology_term_id": "",
    }
    for dim in dims:
        dims[dim] = (
            find_all_dim_option_values(snapshot, criteria.organism_ontology_term_id, dim)
            if is_criteria_empty(criteria)
            else find_dim_option_values(criteria, snapshot, dim)
        )

    # For schema-4 we filter out comma-delimited values for `self_reported_ethnicity_ontology_term_id`
    # from the options list per functional requirements:
    # See: https://github.com/chanzuckerberg/single-cell/issues/596
    ethnicity_term_ids = dims["self_reported_ethnicity_ontology_term_id"]
    dims["self_reported_ethnicity_ontology_term_id"] = [term_id for term_id in ethnicity_term_ids if "," not in term_id]

    response_filter_dims_values = dict(
        disease_terms=build_ontology_term_id_label_mapping(dims["disease_ontology_term_id"]),
        sex_terms=build_ontology_term_id_label_mapping(dims["sex_ontology_term_id"]),
        development_stage_terms=build_ontology_term_id_label_mapping(dims["development_stage_ontology_term_id"]),
        self_reported_ethnicity_terms=build_ontology_term_id_label_mapping(
            dims["self_reported_ethnicity_ontology_term_id"]
        ),
        publication_citations=dims["publication_citation"],
        cell_type_terms=build_ontology_term_id_label_mapping(dims["cell_type_ontology_term_id"]),
        tissue_terms=build_ontology_term_id_label_mapping(dims["tissue_ontology_term_id"]),
        organism_terms=build_ontology_term_id_label_mapping(dims["organism_ontology_term_id"]),
    )

    return response_filter_dims_values


@tracer.wrap(name="build_gene_id_label_mapping", service="wmg-api", resource="query", span_type="wmg-api")
def build_gene_id_label_mapping(gene_ontology_term_ids: List[str]) -> List[dict]:
    return [
        {gene_ontology_term_id: gene_term_label(gene_ontology_term_id)}
        for gene_ontology_term_id in gene_ontology_term_ids
    ]


def build_ontology_term_id_label_mapping(ontology_term_ids: Iterable[str]) -> List[dict]:
    return [{ontology_term_id: ontology_term_label(ontology_term_id)} for ontology_term_id in ontology_term_ids]


@tracer.wrap(name="differentialExpression", service="de-api", resource="differentialExpression", span_type="de-api")
def differentialExpression():
    request = connexion.request.json

    queryGroup1Filters = request["queryGroup1Filters"]
    queryGroup2Filters = request["queryGroup2Filters"]

    criteria1 = DeQueryCriteria(**queryGroup1Filters)
    criteria2 = DeQueryCriteria(**queryGroup2Filters)

    snapshot: WmgSnapshot = load_snapshot(
        snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION,
        explicit_snapshot_id_to_load=WMG_API_FORCE_LOAD_SNAPSHOT_ID,
        snapshot_fs_root_path=SNAPSHOT_FS_ROOT_PATH,
    )

    # cube_query_params are not required to instantiate WmgQuery for differential expression
    q = WmgQuery(snapshot, cube_query_params=None)

    with ServerTiming.time("run differential expression"):
        de_results, query_group_cell_counts, successCode = run_differential_expression(q, criteria1, criteria2)

    return jsonify(
        dict(
            snapshot_id=snapshot.snapshot_identifier,
            differentialExpressionResults=de_results,
            queryGroupCellCounts=query_group_cell_counts,
            successCode=successCode,
        )
    )


def run_differential_expression(q: WmgQuery, criteria1, criteria2) -> Tuple[List[Dict], int]:
    """
    Runs differential expression analysis between two sets of criteria.

    This function takes two criteria objects, representing two different groups for comparison
    in a differential expression analysis. It first augments these criteria with their descendant
    cell types if cell_type_ontology_term_ids are specified. Then, it calculates cell counts for each
    group, filters out overlapping populations, aggregates expression summaries by gene ontology term IDs,
    and finally runs the statistical test (t-test).

    Parameters:
    - q: WmgQuery object
    - criteria1: The first set of criteria for differential expression analysis.
    - criteria2: The second set of criteria for differential expression analysis.

    Returns:
    A tuple containing two elements:
    - A list of dictionaries, each representing a gene and its differential expression metrics.
    - An integer success code, where 1 indicates an issue (e.g., no data after filtering overlapping populations).
    """

    # augment criteria1 and criteria2 with descendants if cell_type_ontology_term_ids is specified
    # this is effectively rollup
    if criteria1.cell_type_ontology_term_ids:
        criteria1.cell_type_ontology_term_ids = list(
            set(sum([descendants(i) for i in criteria1.cell_type_ontology_term_ids], []))
        )
    if criteria2.cell_type_ontology_term_ids:
        criteria2.cell_type_ontology_term_ids = list(
            set(sum([descendants(i) for i in criteria2.cell_type_ontology_term_ids], []))
        )

    cell_counts1 = q.cell_counts(criteria1)
    cell_counts2 = q.cell_counts(criteria2)

    n_cells1 = cell_counts1["n_total_cells"].sum()
    n_cells2 = cell_counts2["n_total_cells"].sum()
    es1 = q.expression_summary_diffexp(criteria1)
    es2 = q.expression_summary_diffexp(criteria2)

    # Skip this step for now, validation may identify that it is required
    # # filter out rows from es2 that are in es1
    # # this prevents overlapping populations from being compared
    # filter_columns = [
    #     col
    #     for col in (base_expression_summary_indexed_dims + expression_summary_secondary_dims)
    #     if col in es1.columns and col in es2.columns
    # ]
    # index1 = es1.set_index(filter_columns).index
    # index2 = es2.set_index(filter_columns).index
    # es2 = es2[~index2.isin(index1)]
    # if es2.shape[0] == 0:
    #     return [], 1

    es_agg1 = es1.groupby("gene_ontology_term_id").sum(numeric_only=True)
    es_agg2 = es2.groupby("gene_ontology_term_id").sum(numeric_only=True)

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

    pvals, effects, tscores = _run_ttest(sums1, sqsums1, n_cells1, sums2, sqsums2, n_cells2)
    de_genes = np.array(genes)[np.argsort(-effects)]
    p = pvals[np.argsort(-effects)]
    tscores = tscores[np.argsort(-effects)]
    effects = effects[np.argsort(-effects)]

    statistics = []
    for i in range(len(p)):
        pi = p[i]
        ei = effects[i]
        ti = tscores[i]
        if ei is not np.nan and pi is not np.nan:
            statistics.append(
                {
                    "gene_ontology_term_id": de_genes[i],
                    "gene_symbol": gene_term_label(de_genes[i]),
                    "p_value": pi,
                    "effect_size": ei,
                    "t_score": ti,
                }
            )

    query_group_cell_counts = {
        "queryGroup1CellCount": int(n_cells1),
        "queryGroup2CellCount": int(n_cells2),
    }
    return statistics, query_group_cell_counts, 0


def _run_ttest(sum1, sumsq1, n1, sum2, sumsq2, n2):
    """
    Calculate t-test statistics.

    Arguments
    ---------
    sum1 - np.ndarray (1 x M)
        Array of sum expression in target pop for each gene
    sumsq1 - np.ndarray (1 x M)
        Array of sum of squared expressions in target pop for each gene
    n1 - np.ndarray (1 x M)
        Number of cells in target pop for each gene
    sum2 - np.ndarray (N x M)
        Array of sum expression in context pops for each gene
    sumsq2 - np.ndarray (N x M)
        Array of sum of squared expressions in context pops for each gene
    n2 - np.ndarray (N x M)
        Number of cells in context pops for each gene

    Returns
    -------
    pvals_adj - np.ndarray
        adjusted p-values for each comparison for each gene
    effects - np.ndarray
        effect sizes for each comparison for each gene
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        mean1 = sum1 / n1
        meansq1 = sumsq1 / n1

        mean2 = sum2 / n2
        meansq2 = sumsq2 / n2

        var1 = meansq1 - mean1**2
        var1[var1 < 0] = 0
        var2 = meansq2 - mean2**2
        var2[var2 < 0] = 0

        var1_n = var1 / n1
        var2_n = var2 / n2
        sum_var_n = var1_n + var2_n
        dof = sum_var_n**2 / (var1_n**2 / (n1 - 1) + var2_n**2 / (n2 - 1))
        tscores = (mean1 - mean2) / np.sqrt(sum_var_n)
        effects = (mean1 - mean2) / np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 1))

    pvals = stats.t.sf(tscores, dof)
    return pvals, effects, tscores
