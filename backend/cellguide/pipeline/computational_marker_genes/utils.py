import os
from typing import Tuple

import numpy as np
from numba import njit, prange

from backend.common.constants import DEPLOYMENT_STAGE_TO_API_URL
from backend.wmg.data.utils import setup_retry_session


@njit(parallel=True)
def calculate_specificity_excluding_nans(treatment, control):
    treatment = treatment.flatten()

    specificities = np.zeros(treatment.size)
    for i in prange(treatment.size):
        if np.isnan(treatment[i]):
            continue
        col = control[:, i]
        col = col[~np.isnan(col)]
        specificities[i] = (treatment[i] > col).mean()
    return specificities


def calculate_cohens_d(
    *, sum1: np.ndarray, sumsq1: np.ndarray, n1: np.ndarray, sum2: np.ndarray, sumsq2: np.ndarray, n2: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculates Cohen's d for two sets of data.
    Arrays "1" and "2" have to be broadcastable into each other.

    Arguments
    ---------
    sum1 - Sum of the first set of data
    sumsq1 - Sum of squares of the first set of data
    n1 - Number of elements in the first set of data
    sum2 - Sum of the second set of data
    sumsq2 - Sum of squares of the second set of data
    n2 - Number of elements in the second set of data

    Returns
    -------
    pvals - The p-values of the t-test
    effects - The effect sizes of the t-test
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
        effects = (mean1 - mean2) / np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 1))

    return effects


def query_gene_info_for_gene_description(gene_id: str) -> str:
    """
    Query the gene description for a given gene ID.

    Arguments
    ---------
    gene_id - The ID of the gene to query

    Returns
    -------
    The name of the gene if the query is successful, otherwise returns the gene ID.
    """

    deployment_stage = os.environ.get("DEPLOYMENT_STAGE")
    # if deployment stage is rdev, use staging API for this query
    # otherwise, we run into JSON decode errors because the gene info endpoint
    # is behind Auth0.
    API_URL = DEPLOYMENT_STAGE_TO_API_URL.get(
        deployment_stage, "https://api.cellxgene.staging.single-cell.czi.technology"
    )

    session = setup_retry_session()
    response = session.get(f"{API_URL}/gene_info/v1/gene_info?gene={gene_id}")
    if response.status_code == 200:
        data = response.json()
        return data["name"]
    else:
        return gene_id


@njit(parallel=True)
def bootstrap_rows_percentiles(
    X: np.ndarray, num_replicates: int = 1000, num_samples: int = 100, percentile: float = 5
):
    """
    This function bootstraps rows of a given matrix X.

    Arguments
    ---------
    X : np.ndarray
        The input matrix to bootstrap.
    num_replicates : int, optional
        The number of bootstrap replicates to generate, by default 1000.
    num_samples : int, optional
        The number of samples to draw in each bootstrap replicate, by default 100.
    percentile : float, optional
        The percentile of the bootstrapped samples for each replicate, by default 15.

    Returns
    -------
    bootstrap_percentile : np.ndarray
        The percentile of the bootstrapped samples for each replicate.
    """

    bootstrap_percentile = np.zeros((num_replicates, X.shape[1]), dtype="float")

    indices = np.random.randint(0, X.shape[0], size=(num_replicates, num_samples))
    # for each replicate
    for n_i in prange(num_replicates):
        bootstrap_percentile[n_i] = sort_matrix_columns(X[indices[n_i]], percentile, num_samples)

    return bootstrap_percentile


@njit
def sort_matrix_columns(matrix, percentile, num_samples):
    """
    This function sorts the columns of a given matrix and returns the index associated with
    the specified percentile of the sorted samples for each column. This approximates
    np.nanpercentile(matrix, percentile, axis=0).

    Arguments
    ---------
    matrix : np.ndarray
        The input matrix to sort.
    percentile : float
        The percentile of the sorted samples for each column.
    num_samples : int
        The number of samples in each column.

    Returns
    -------
    result : np.ndarray
        The sorted columns of the input matrix.
    """
    num_cols = matrix.shape[1]
    result = np.empty(num_cols)
    for col in range(num_cols):
        sorted_col = np.sort(matrix[:, col])
        num_nans = np.isnan(sorted_col).sum()
        num_non_nans = num_samples - num_nans
        sample_index = int(np.round(percentile / 100 * num_non_nans))
        result[col] = sorted_col[sample_index]
    return result
