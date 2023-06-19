import json
import logging
import os
import sys
from typing import Dict, List

import cellxgene_schema

from backend.layers.business.business import BusinessLogic
from backend.layers.business.entities import CollectionQueryFilter
from backend.layers.common.entities import DatasetArtifactType
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.processing import logger
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.uri_provider import UriProvider

logger.configure_logging(level=logging.INFO)


database_provider = DatabaseProvider()
s3_provider = S3Provider()
uri_provider = UriProvider()

business_logic = BusinessLogic(
    database_provider,
    None,  # Not required
    None,  # Not required
    s3_provider,
    uri_provider,
)


def gather_collections() -> Dict[str, Dict[str, List[str]]]:
    """
    This function is used to gather all the collections and their datasets that will be migrated
    :return: A dictionary with the following structure:
    {
        "published": {"collection_id": ["dataset_id", "dataset_id", ...], ...},
        "revision": {"collection_version_id": ["dataset_id", "dataset_id", ...], ...}
        "private": {"collection_id": ["dataset_id", "dataset_id", ...], ...},
    }
    """

    collections_filter = CollectionQueryFilter()
    collections = business_logic.get_collections(collections_filter)

    response = {"published": {}, "private": {}, "revision": {}}

    def get_datasets(_version):
        datasets = []
        for dataset in _version.datasets:
            datasets.append({"dataset_id": dataset.dataset_id.id, "dataset_version_id": dataset.version_id.id})
        return datasets

    has_revision = []  # list of collections to skip if published with an active revision
    for collection in collections:
        if collection.is_published() and collection.collection_id not in has_revision:
            version = business_logic.get_collection_version(collection.version_id)
            response["published"][version.collection_id.id] = get_datasets(version)
        elif collection.is_unpublished_version():
            has_revision.append(
                collection.collection_id
            )  # skip the published version if it hasn't been encountered yet
            response["published"].pop(collection.collection_id.id, None)
            # remove from published if it was already encountered
            version = business_logic.get_collection_version(collection.version_id)
            response["revision"][version.version_id.id] = get_datasets(version)
            # using version id instead of collection id
        elif collection.is_initial_unpublished_version():
            version = business_logic.get_collection_version(collection.version_id)
            response["private"][version.collection_id.id] = get_datasets(version)
    return response


def migrate(step_name: str) -> bool:
    """
    Gets called by the step function at every different step, as defined by `step_name`
    """
    if step_name == "corpus_migrate":
        corpus_migrate()
    if step_name == "collection_migrate":
        collection_id = os.environ["collection_id"]
        datasets = json.loads(os.environ["datasets"])
        can_open_revision = bool(os.environ["can_open_revision"])  # TODO: Check syntax for str -> bool
        collection_migrate(collection_id, datasets, can_open_revision)
    if step_name == "dataset_migrate":
        collection_id = os.environ["collection_id"]
        dataset_id = os.environ["dataset_id"]
        dataset_version_id = os.envion["dataset_version_id"]
        dataset_migrate(collection_id, dataset_id, dataset_version_id)
    return True


def dataset_migrate(collection_id: str, dataset_id: str, dataset_version_id: str) -> Dict[str, str]:
    raw_h5ad_uri = [
        artifact.uri
        for artifact in business_logic.get_dataset_artifacts(dataset_version_id)
        if artifact.type == DatasetArtifactType.RAW_H5AD
    ]
    bucket_name, object_key = s3_provider.parse_s3_uri(raw_h5ad_uri)
    s3_provider.download_file(bucket_name, object_key, "raw.h5ad")
    cellxgene_schema.migrate("raw.h5ad", "migrated.h5ad", collection_id, dataset_id)
    # TODO: define where this upload bucket actually should be
    upload_bucket_name = os.environ["UPLOAD_BUCKET"]
    dst_uri = f"{dataset_version_id}/migrated.h5ad"
    s3_provider.upload_file("migrated.h5ad", upload_bucket_name, dst_uri)
    url = f"s3://{upload_bucket_name}/{dst_uri}"
    new_dataset_version_id, _ = business_logic.ingest_dataset(collection_id, url, None, dataset_version_id)
    return {"new_dataset_version_id": new_dataset_version_id}


def collection_migrate(
    collection_id: str, datasets: List[Dict[str, str]], can_open_revision: bool
) -> List[Dict[str, str]]:
    private_collection_id = collection_id
    if can_open_revision:
        private_collection_id = business_logic.create_collection_version(collection_id).version_id.id
    return [
        {
            "collection_id": private_collection_id,
            "dataset_id": dataset["dataset_id"],
            "dataset_version_id": dataset["dataset_version_id"],
        }
        for dataset in datasets
    ]


def corpus_migrate() -> bool:
    gather_collections()
    return True


if __name__ == "__main__":
    step_name = os.environ["STEP_NAME"]
    rv = migrate(step_name)
    sys.exit(rv)
