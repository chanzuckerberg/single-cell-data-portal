import os
from backend.layers.common.entities import DatasetVersion


def generate(dataset: DatasetVersion, use_canonical=True):
    """
    Generates an explorer_url for the present dataset.
    By default, set to use canonical dataset ID to build URL but
    will use dataset version ID if use_canonical is set to False.
    """
    frontend_url = os.getenv("FRONTEND_URL", "")
    dataset_id = dataset.dataset_id if use_canonical else dataset.version_id
    return f"{frontend_url}/e/{dataset_id}.cxg/"
