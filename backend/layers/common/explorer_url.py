import os
from backend.layers.common.entities import DatasetVersion

def generate(dataset: DatasetVersion, is_published=False):
    """
    Generates an explorer_url for the present dataset.
    If the dataset is published, or if it's referenced in the context of a published collection,
    the dataset_id should be used. Otherwise, version_id should be used.
    """
    frontend_url = os.getenv("FRONTEND_URL")
    id = dataset.dataset_id if is_published else dataset.version_id
    return f"{frontend_url}/{id}.cxg"