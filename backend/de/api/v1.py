import os
from typing import Any, Dict, Iterable, List, Tuple

import connexion
import numpy as np
import pandas as pd
from ddtrace import tracer
from flask import jsonify
from scipy import stats
from server_timing import Timing as ServerTiming

from backend.common.census_cube.data.criteria import BaseQueryCriteria
from backend.common.census_cube.data.ontology_labels import gene_term_label, ontology_term_label
from backend.common.census_cube.data.query import CensusCubeQuery
from backend.common.census_cube.data.schemas.cube_schema_diffexp import cell_counts_logical_dims_exclude_dataset_id
from backend.common.census_cube.data.snapshot import CensusSnapshot, load_snapshot
from backend.common.marker_gene_files.blacklist import marker_gene_blacklist
from backend.common.utils.rollup import descendants
from backend.wmg.api.wmg_api_config import (
    WMG_API_FORCE_LOAD_SNAPSHOT_ID,
    WMG_API_READ_FS_CACHED_SNAPSHOT,
    WMG_API_SNAPSHOT_FS_CACHE_ROOT_PATH,
    WMG_API_SNAPSHOT_SCHEMA_VERSION,
)

DEPLOYMENT_STAGE = os.environ.get("DEPLOYMENT_STAGE", "")
SNAPSHOT_FS_ROOT_PATH = (
    WMG_API_SNAPSHOT_FS_CACHE_ROOT_PATH if (WMG_API_READ_FS_CACHED_SNAPSHOT and DEPLOYMENT_STAGE != "test") else None
)

"""
We have a data structure that must be shared between CellGuide, DE, and WMG applications as they all are derived from the same source.

Which means that any WMG data stuff must be generalized. What's a better name for it? WMG stands for "Where's my Gene" but it's vestigial, 
the new name is Gene Expression. The data structure contains expression aggregates across cells, but now it also contains other information to serve
cases like marker genes and differential expression.

A better name would be: 
"""


@tracer.wrap(name="filters", service="wmg-api", resource="filters", span_type="wmg-api")
def filters():
    request = connexion.request.json
    sanitize_api_query_dict(request["filter"])

    criteria = BaseQueryCriteria(**request["filter"])

    with ServerTiming.time("load snapshot"):
        snapshot: CensusSnapshot = load_snapshot(
            snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION,
            explicit_snapshot_id_to_load=WMG_API_FORCE_LOAD_SNAPSHOT_ID,
            snapshot_fs_root_path=SNAPSHOT_FS_ROOT_PATH,
        )

    with ServerTiming.time("calculate filters and build response"):
        q = CensusCubeQuery(snapshot, cube_query_params=None)
        response_filter_dims_values = build_filter_dims_values(criteria, snapshot, q)
        n_cells = _get_cell_counts_for_query(q, criteria)

        response = jsonify(
            dict(
                snapshot_id=snapshot.snapshot_identifier,
                filter_dims=response_filter_dims_values,
                n_cells=n_cells,
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


def is_criteria_empty(criteria: BaseQueryCriteria) -> bool:
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
def build_filter_dims_values(criteria: BaseQueryCriteria, snapshot: CensusSnapshot, q: CensusCubeQuery) -> Dict:

    if is_criteria_empty(criteria):
        df = snapshot.cell_counts_df[
            snapshot.cell_counts_df["organism_ontology_term_id"] == criteria.organism_ontology_term_id
        ]
        dims = {col: df[col].unique().tolist() for col in df.columns}
    else:
        df_ref = q.cell_counts_df(criteria)
        dims = {}
        for key in criteria.dict():
            col_name = key[:-1] if key[-1] == "s" else key
            if key == "organism_ontology_term_id" or col_name not in df_ref:
                continue

            values = getattr(criteria, key)
            if len(values) == 0:
                dims[col_name] = df_ref[col_name].unique().tolist()
            else:
                setattr(criteria, key, [])

                df = q.cell_counts_df(criteria)
                dims[col_name] = df[col_name].unique().tolist()

                setattr(criteria, key, values)

        dims["organism_ontology_term_id"] = snapshot.cell_counts_df["organism_ontology_term_id"].unique().tolist()

    # For schema-4 we filter out comma-delimited values for `self_reported_ethnicity_ontology_term_id`
    # from the options list per functional requirements:
    # See: https://github.com/chanzuckerberg/single-cell/issues/596
    ethnicity_term_ids = dims["self_reported_ethnicity_ontology_term_id"]
    dims["self_reported_ethnicity_ontology_term_id"] = [term_id for term_id in ethnicity_term_ids if "," not in term_id]

    response_filter_dims_values = dict(
        datasets=fetch_datasets_metadata(snapshot, dims["dataset_id"]),
        disease_terms=build_ontology_term_id_label_mapping(dims["disease_ontology_term_id"]),
        sex_terms=build_ontology_term_id_label_mapping(dims["sex_ontology_term_id"]),
        development_stage_terms=[],
        self_reported_ethnicity_terms=build_ontology_term_id_label_mapping(
            dims["self_reported_ethnicity_ontology_term_id"]
        ),
        publication_citations=dims["publication_citation"],
        cell_type_terms=build_ontology_term_id_label_mapping(dims["cell_type_ontology_term_id"]),
        tissue_terms=build_ontology_term_id_label_mapping(dims["tissue_ontology_term_id"]),
        organism_terms=build_ontology_term_id_label_mapping(dims["organism_ontology_term_id"]),
    )

    return response_filter_dims_values


def fetch_datasets_metadata(snapshot: CensusSnapshot, dataset_ids: Iterable[str]) -> List[Dict]:
    return [
        snapshot.dataset_metadata.get(dataset_id, dict(id=dataset_id, label="", collection_id="", collection_label=""))
        for dataset_id in dataset_ids
    ]


@tracer.wrap(name="build_gene_id_label_mapping", service="wmg-api", resource="query", span_type="wmg-api")
def build_gene_id_label_mapping(gene_ontology_term_ids: List[str]) -> List[dict]:
    return [
        {gene_ontology_term_id: gene_term_label(gene_ontology_term_id)}
        for gene_ontology_term_id in gene_ontology_term_ids
    ]


def build_ontology_term_id_label_mapping(ontology_term_ids: Iterable[str]) -> List[dict]:
    return [ontology_term_id_label_mapping(ontology_term_id) for ontology_term_id in ontology_term_ids]


def ontology_term_id_label_mapping(ontology_term_id: str) -> dict:
    return {ontology_term_id: ontology_term_label(ontology_term_id)}


@tracer.wrap(name="differentialExpression", service="de-api", resource="differentialExpression", span_type="de-api")
def differentialExpression():
    request = connexion.request.json

    queryGroup1Filters = request["queryGroup1Filters"]
    queryGroup2Filters = request["queryGroup2Filters"]

    criteria1 = BaseQueryCriteria(**queryGroup1Filters)
    criteria2 = BaseQueryCriteria(**queryGroup2Filters)

    snapshot: CensusSnapshot = load_snapshot(
        snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION,
        explicit_snapshot_id_to_load=WMG_API_FORCE_LOAD_SNAPSHOT_ID,
        snapshot_fs_root_path=SNAPSHOT_FS_ROOT_PATH,
    )

    # cube_query_params are not required to instantiate CensusCubeQuery for differential expression
    q = CensusCubeQuery(snapshot, cube_query_params=None)

    with ServerTiming.time("run differential expression"):
        de_results, n_overlap = run_differential_expression(q, criteria1, criteria2)

    return jsonify(
        dict(
            snapshot_id=snapshot.snapshot_identifier,
            differentialExpressionResults=de_results,
            n_overlap=n_overlap,
        )
    )


def run_differential_expression(q: CensusCubeQuery, criteria1, criteria2) -> Tuple[List[Dict], int]:
    """
    Runs differential expression analysis between two sets of criteria.

    This function takes two criteria objects, representing two different groups for comparison
    in a differential expression analysis. It first augments these criteria with their descendant
    cell types if cell_type_ontology_term_ids are specified. Then, it calculates cell counts for each
    group, filters out overlapping populations, aggregates expression summaries by gene ontology term IDs,
    and finally runs the statistical test (t-test).

    Parameters:
    - q: CensusCubeQuery object
    - criteria1: The first set of criteria for differential expression analysis.
    - criteria2: The second set of criteria for differential expression analysis.

    Returns:
    A tuple containing two elements:
    - A list of dictionaries, each representing a gene and its differential expression metrics.
    - An integer representing the number of overlapping populations between the two groups.
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

    es1, cell_counts1 = q.expression_summary_and_cell_counts_diffexp(criteria1)
    es2, cell_counts2 = q.expression_summary_and_cell_counts_diffexp(criteria2)

    n_cells1 = cell_counts1["n_total_cells"].sum()
    n_cells2 = cell_counts2["n_total_cells"].sum()

    # identify number of overlapping populations
    filter_columns = [
        col
        for col in cell_counts_logical_dims_exclude_dataset_id
        if col in cell_counts1.columns and col in cell_counts2.columns
    ]
    index1 = cell_counts1.set_index(filter_columns).index
    index2 = cell_counts2.set_index(filter_columns).index
    overlap_filter = index1.isin(index2)
    n_overlap = int(cell_counts1[overlap_filter]["n_total_cells"].sum())

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

    lfc, effects, pvals_adj = _calculate_t_test_metrics(sums1, sqsums1, n_cells1, sums2, sqsums2, n_cells2)
    de_genes = np.array(genes)[np.argsort(-effects)]
    lfc = lfc[np.argsort(-effects)]
    pvals_adj = pvals_adj[np.argsort(-effects)]
    effects = effects[np.argsort(-effects)]

    statistics = []
    for i in range(len(lfc)):
        ei = effects[i]
        pval = pvals_adj[i]
        if ei is not np.nan and pval is not np.nan and de_genes[i] not in marker_gene_blacklist:
            statistics.append(
                {
                    "gene_ontology_term_id": de_genes[i],
                    "gene_symbol": gene_term_label(de_genes[i]),
                    "effect_size": ei,
                    "log_fold_change": lfc[i],
                    "adjusted_p_value": pval,
                }
            )
    return statistics, n_overlap


def _get_cell_counts_for_query(q: CensusCubeQuery, criteria: BaseQueryCriteria) -> pd.DataFrame:
    if criteria.cell_type_ontology_term_ids:
        criteria.cell_type_ontology_term_ids = list(
            set(sum([descendants(i) for i in criteria.cell_type_ontology_term_ids], []))
        )
    cell_counts = q.cell_counts_diffexp_df(criteria)
    return int(cell_counts["n_total_cells"].sum())


def _calculate_t_test_metrics(sum1, sumsq1, n1, sum2, sumsq2, n2):
    """
    Calculate log fold changes, effect sizes, and BH-adjusted p-values for each gene.

    Arguments
    ---------
    sum1 - np.ndarray (1 x M)
        Array of sum expression for each gene in pop 1
    sumsq1 - np.ndarray (1 x M)
        Array of sum of squared expressions for each gene in pop 1
    n1 - int
        Number of cells in target pop for each gene
    sum2 - np.ndarray (1 x M)
        Array of sum expression for each gene in pop 2
    sumsq2 - np.ndarray (1 x M)
        Array of sum of squared expressions for each gene in pop 2
    n2 - int
        Number of cells in pop 2 for each gene

    Returns
    -------
    log_fold_changes - np.ndarray
        log fold changes for each gene
    effects - np.ndarray
        effect sizes for each gene
    pvals_adj - np.ndarray
        adjusted p-values for each gene
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

        # this is log fold change because mean1 and mean2 are already log transformed
        log_fold_changes = mean1 - mean2

        # two-sided test
        pvals = 2 * stats.t.sf(np.abs(tscores), dof)
        pvals_adj = _benjamini_hochberg_correction(pvals)

    return log_fold_changes, effects, pvals_adj


def _benjamini_hochberg_correction(pvals):
    """
    Perform Benjamini-Hochberg correction for multiple testing in a vectorized manner.

    Arguments
    ---------
    pvals - np.ndarray
        Array of p-values to correct

    Returns
    -------
    pvals_adj - np.ndarray
        Adjusted p-values
    """
    n = len(pvals)
    sorted_indices = np.argsort(pvals)
    sorted_pvals = pvals[sorted_indices]
    adjusted_pvals = np.zeros(n)
    cumulative_min = np.minimum.accumulate((sorted_pvals * n / (np.arange(n) + 1))[::-1])[::-1]
    adjusted_pvals[sorted_indices] = cumulative_min
    return adjusted_pvals
