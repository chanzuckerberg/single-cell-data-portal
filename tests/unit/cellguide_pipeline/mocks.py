import json

import numpy as np

from tests.unit.cellguide_pipeline.constants import (
    ASCTB_MASTER_SHEET_FIXTURE_FILENAME,
    CELLGUIDE_PIPELINE_FIXTURES_BASEPATH,
    FAKE_COLLECTIONS_FIXTURE_FILENAME,
    FAKE_DATASETS_FIXTURE_FILENAME,
)


def mock_get_asctb_master_sheet():
    with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{ASCTB_MASTER_SHEET_FIXTURE_FILENAME}", "r") as f:
        return json.load(f)


def mock_get_datasets_from_curation_endpoint():
    with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{FAKE_DATASETS_FIXTURE_FILENAME}", "r") as f:
        return json.load(f)


def mock_get_collections_from_curation_endpoint():
    with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{FAKE_COLLECTIONS_FIXTURE_FILENAME}", "r") as f:
        return json.load(f)


def mock_get_title_and_citation_from_doi(doi: str):
    return f"title from {doi}\n\n-- citation"


def mock_bootstrap_rows_percentiles(
    X: np.ndarray,
    random_indices: np.ndarray,
    num_replicates: int = 1000,
    num_samples: int = 100,
    percentile: float = 5,
):
    """
    Mock the bootstrapping function to return deterministic results.

    (alec) I attempted seeding the RNG but that did not guarantee deterministic results
    for whatever reason. This is a hacky workaround to unblock the unit test.
    """
    return np.tile(np.nanpercentile(X, percentile, axis=0)[None, :], (num_replicates, 1))
