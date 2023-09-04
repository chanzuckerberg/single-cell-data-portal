import json
from typing import Dict, Iterable, List

import connexion
import numpy as np
import pandas as pd
from flask import jsonify
from scipy import stats
from scipy.stats import rankdata
from server_timing import Timing as ServerTiming

from backend.de.api.openai_utils import get_query_from_user_input
from backend.de.data.ontology_labels import gene_term_label, ontology_term_label
from backend.de.data.query import (
    DeQuery,
    DeQueryCriteria,
)
from backend.de.data.snapshot import DeSnapshot, load_snapshot
from backend.de.data.utils import find_all_dim_option_values, find_dim_option_values

# load gene sets for gene set enrichment analysis


def _load_gene_sets_for_enrichment_analysis():
    with open("backend/de/data/fixtures/REACTOME_2022_GENE_SETS.json") as f:
        global GENE_SETS_FOR_ENRICHMENT_ANALYSIS
        GENE_SETS_FOR_ENRICHMENT_ANALYSIS = json.load(f)


_load_gene_sets_for_enrichment_analysis()


def filters():
    request = connexion.request.json
    criteria = DeQueryCriteria(**request["filter"])

    with ServerTiming.time("calculate filters and build response"):
        snapshot: DeSnapshot = load_snapshot()
        response_filter_dims_values = build_filter_dims_values(criteria, snapshot)
        response = jsonify(
            dict(
                snapshot_id=snapshot.snapshot_identifier,
                filter_dims=response_filter_dims_values,
            )
        )
    return response


REMAP_DIMENSION_NAMES = {
    "disease_ontology_term_ids": "disease_terms",
    "cell_type_ontology_term_ids": "cell_type_terms",
    "tissue_ontology_term_ids": "tissue_terms",
    "organism_ontology_term_id": "organism_terms",
    "development_stage_ontology_term_ids": "development_stage_terms",
    "sex_ontology_term_ids": "sex_terms",
    "self_reported_ethnicity_ontology_term_ids": "self_reported_ethnicity_terms",
}


def remap_dimension_names(criteria: dict[str, list[str]]) -> dict[str, list[str]]:
    return {REMAP_DIMENSION_NAMES.get(k, k): build_ontology_term_id_label_mapping(v) for k, v in criteria.items()}


def getDeQuery():
    request = connexion.request.json
    user_query = request["user_query"]

    with ServerTiming.time("translate user input to query criteria"):
        snapshot: DeSnapshot = load_snapshot()
        query_criteria1, query_criteria2 = get_query_from_user_input(user_query, snapshot)

        response = jsonify(
            dict(
                snapshot_id=snapshot.snapshot_identifier,
                query_criteria1=remap_dimension_names(query_criteria1),
                query_criteria2=remap_dimension_names(query_criteria2),
            )
        )
    return response


def differentialExpression():
    global GENE_SETS_FOR_ENRICHMENT_ANALYSIS

    request = connexion.request.json

    queryGroup1Filters = request["queryGroup1Filters"]
    queryGroup2Filters = request["queryGroup2Filters"]

    criteria1 = DeQueryCriteria(**queryGroup1Filters)
    criteria2 = DeQueryCriteria(**queryGroup2Filters)

    snapshot: DeSnapshot = load_snapshot()

    q = DeQuery(snapshot)

    with ServerTiming.time("run differential expression"):
        results1, results2 = run_differential_expression(q, criteria1, criteria2)

    with ServerTiming.time("run enrichment analysis"):
        ordered_genes1 = list(pd.DataFrame(results1).sort_values("t_score", ascending=False)["gene_symbol"])
        ordered_genes2 = list(pd.DataFrame(results2).sort_values("t_score", ascending=False)["gene_symbol"])
        goea_results1 = GOEA(ordered_genes1, GENE_SETS_FOR_ENRICHMENT_ANALYSIS)
        goea_results2 = GOEA(ordered_genes2, GENE_SETS_FOR_ENRICHMENT_ANALYSIS)

    return jsonify(
        dict(
            snapshot_id=snapshot.snapshot_identifier,
            differentialExpressionResults1=results1,
            differentialExpressionResults2=results2,
            pathwayEnrichmentAnalysisResults1=goea_results1.to_dict(orient="records"),
            pathwayEnrichmentAnalysisResults2=goea_results2.to_dict(orient="records"),
            successCode=0,
        )
    )


def fetch_datasets_metadata(snapshot: DeSnapshot, dataset_ids: Iterable[str]) -> List[Dict]:
    return [
        snapshot.dataset_dict.get(dataset_id, dict(id=dataset_id, label="", collection_id="", collection_label=""))
        for dataset_id in dataset_ids
    ]


def is_criteria_empty(criteria: DeQueryCriteria) -> bool:
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


def build_filter_dims_values(criteria: DeQueryCriteria, snapshot: DeSnapshot) -> Dict:
    dims = {
        "dataset_id": "",
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
        organism_terms=build_ontology_term_id_label_mapping(dims["organism_ontology_term_id"]),
    )

    return response_filter_dims_values


def build_ontology_term_id_label_mapping(ontology_term_ids: Iterable[str]) -> List[dict]:
    return [{ontology_term_id: ontology_term_label(ontology_term_id)} for ontology_term_id in ontology_term_ids]


def run_differential_expression(q, criteria1, criteria2, pval_thr=1e-5, n_genes=500):
    cell_counts1 = q.cell_counts(criteria1)
    cell_counts2 = q.cell_counts(criteria2)

    n_cells1 = cell_counts1["n_total_cells"].sum()
    n_cells2 = cell_counts2["n_total_cells"].sum()
    es_agg1 = q.expression_summary(criteria1).groupby("gene_ontology_term_id").sum(numeric_only=True)
    es_agg2 = q.expression_summary(criteria2).groupby("gene_ontology_term_id").sum(numeric_only=True)
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

    statistics1 = []
    for i in range(len(p)):
        pi = p[i]
        ei = effects[i]
        ti = tscores[i]
        if ei is not np.nan and pi is not np.nan and pi < pval_thr:
            statistics1.append(
                {
                    "gene_ontology_term_id": de_genes[i],
                    "gene_symbol": gene_term_label(de_genes[i]),
                    "p_value": pi,
                    "effect_size": ei,
                    "t_score": ti,
                }
            )
            if len(statistics1) >= n_genes and n_genes > 0:
                break

    pvals, effects, tscores = _run_ttest(sums2, sqsums2, n_cells2, sums1, sqsums1, n_cells1)
    de_genes = np.array(genes)[np.argsort(-effects)]
    p = pvals[np.argsort(-effects)]
    tscores = tscores[np.argsort(-effects)]
    effects = effects[np.argsort(-effects)]

    statistics2 = []
    for i in range(len(p)):
        pi = p[i]
        ei = effects[i]
        ti = tscores[i]
        if ei is not np.nan and pi is not np.nan and pi < pval_thr:
            statistics2.append(
                {
                    "gene_ontology_term_id": de_genes[i],
                    "gene_symbol": gene_term_label(de_genes[i]),
                    "p_value": pi,
                    "effect_size": ei,
                    "t_score": ti,
                }
            )
            if len(statistics2) >= n_genes and n_genes > 0:
                break
    return statistics1, statistics2


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


def GOEA(target_genes, GENE_SETS, goterms=None, fdr_thresh=0.25, p_thresh=1e-3):
    target_genes = np.array(target_genes)

    GENE_SETS_IN_TARGET = {i: list(set(GENE_SETS[i]).intersection(target_genes)) for i in GENE_SETS}

    all_genes = np.unique(np.concatenate(list(GENE_SETS.values())))
    all_genes = np.array(all_genes)

    # if goterms is None, use all the goterms found in `GENE_SETS`
    if goterms is None:
        goterms = np.unique(list(GENE_SETS.keys()))
    else:
        goterms = goterms[np.in1d(goterms, np.unique(list(GENE_SETS.keys())))]

    # ensure that target genes are all present in `all_genes`
    _, ix = np.unique(target_genes, return_index=True)
    target_genes = target_genes[np.sort(ix)]
    target_genes = target_genes[np.in1d(target_genes, all_genes)]

    # N -- total number of genes
    N = all_genes.size

    probs = []
    probs_genes = []

    for goterm in goterms:
        # identify genes associated with this go term
        gene_set = GENE_SETS[goterm]
        B = len(gene_set)

        gene_set_in_target = GENE_SETS_IN_TARGET[goterm]
        b = len(gene_set_in_target)

        if b != 0:
            # calculate the enrichment probability as the cumulative sum of the tail end of a hypergeometric distribution
            # with parameters (N,B,n,b)
            n = target_genes.size
            num_iter = min(n, B)
            rng = np.arange(b, num_iter + 1)

            log_factorial = np.arange(max([n, N, B]) + 1).astype("float")
            log_factorial[1:] = np.log(log_factorial[1:]).cumsum()

            def log_binomial(n, k, log_factorial):
                return log_factorial[n] - (log_factorial[k] + log_factorial[n - k])

            probs.append(
                sum(
                    [
                        np.exp(
                            log_binomial(n, i, log_factorial)
                            + log_binomial(N - n, B - i, log_factorial)
                            - log_binomial(N, B, log_factorial)
                        )
                        for i in rng
                    ]
                )
            )
        else:
            probs.append(1.0)

        # append associated genes to a list
        probs_genes.append(gene_set_in_target)

    probs = np.array(probs)
    probs_genes = np.array([";".join(g) for g in probs_genes])

    # adjust p value to correct for multiple testing
    fdr_q_probs = probs.size * probs / rankdata(probs, method="ordinal")

    # filter out go terms based on the FDR q value and p value thresholds
    filt = np.logical_and(fdr_q_probs < fdr_thresh, probs < p_thresh)
    enriched_goterms = goterms[filt]
    p_values = probs[filt]
    fdr_q_probs = fdr_q_probs[filt]
    probs_genes = probs_genes[filt]

    # construct the Pandas DataFrame
    enriched_goterms_df = pd.DataFrame(data=fdr_q_probs, columns=["fdr_q_value"])
    enriched_goterms_df["gene_set"] = enriched_goterms
    enriched_goterms_df["p_value"] = p_values
    enriched_goterms_df["gene_symbols"] = probs_genes

    # sort in ascending order by the p value
    enriched_goterms_df = enriched_goterms_df.sort_values("fdr_q_value")
    return enriched_goterms_df
