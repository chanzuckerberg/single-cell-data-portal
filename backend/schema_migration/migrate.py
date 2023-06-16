import itertools
import logging
import os
import sys
from typing import Dict, List

from backend.layers.business.business import BusinessLogic
from backend.layers.business.entities import CollectionQueryFilter
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

    CollectionQueryFilter()
    published_collections = business_logic.get_collections(CollectionQueryFilter(is_published=True))
    unpublished_collections = business_logic.get_collections(CollectionQueryFilter(is_published=False))

    response = {"published": {}, "private": {}, "revision": {}}

    def get_datasets(_version):
        datasets = []
        for dataset in _version.datasets:
            datasets.append(dataset.dataset_id.id)
        return datasets

    collections = itertools.chain(unpublished_collections, published_collections)
    # evaluate unpublished collections first, so that published versions are skipped if there is an active revision
    has_revision = []  # list of collections to skip if published with an active revision
    for collection in collections:
        if collection.is_published() and collection.collection_id not in has_revision:
            version = business_logic.get_collection_version(collection.version_id)
            response["published"][version.collection_id.id] = get_datasets(version)
        elif collection.is_unpublished_version():
            has_revision.append(collection.collection_id)  # revision found, skip published version
            version = business_logic.get_collection_version(collection.version_id)
            response["revision"][version.version_id.id] = get_datasets(version)
            # using version id instead of collection id
        elif collection.is_initial_unpublished_version():
            version = business_logic.get_collection_version(collection.version_id)
            response["private"][version.collection_id.id] = get_datasets(version)
    return response


def migrate(step_name: str, business_logic: BusinessLogic) -> bool:
    """
    Gets called by the step function at every different step, as defined by `step_name`
    """
    if step_name == "corpus_migrate":
        corpus_migrate(business_logic)
    return True


def corpus_migrate(business_logic: BusinessLogic) -> bool:
    gather_collections()
    return True


def main():
    step_name = os.environ["STEP_NAME"]
    rv = migrate(step_name, business_logic)
    sys.exit(rv)


if __name__ == "__main__":
    main()
