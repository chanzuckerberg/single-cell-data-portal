import itertools
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


class SchemaMigrate:
    def __init__(self, business_logic: BusinessLogic, step_name: str):
        self.business_logic = business_logic
        self.step_name = step_name

    def gather_collections(self) -> Dict[str, Dict[str, List[str]]]:
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
        published_collections = self.business_logic.get_collections(CollectionQueryFilter(is_published=True))
        unpublished_collections = self.business_logic.get_collections(CollectionQueryFilter(is_published=False))

        response = {"published": {}, "private": {}, "revision": {}}

        def get_datasets(_version):
            datasets = []
            for dataset in _version.datasets:
                datasets.append({"dataset_id": dataset.dataset_id.id, "dataset_version_id": dataset.version_id.id})
            return datasets

        collections = itertools.chain(unpublished_collections, published_collections)
        # evaluate unpublished collections first, so that published versions are skipped if there is an active revision
        has_revision = []  # list of collections to skip if published with an active revision
        for collection in collections:
            if collection.is_published() and collection.collection_id not in has_revision:
                version = self.business_logic.get_collection_version(collection.version_id)
                response["published"][version.collection_id.id] = get_datasets(version)
            elif collection.is_unpublished_version():
                has_revision.append(collection.collection_id)  # revision found, skip published version
                version = self.business_logic.get_collection_version(collection.version_id)
                response["revision"][version.version_id.id] = get_datasets(version)
                # using version id instead of collection id
            elif collection.is_initial_unpublished_version():
                version = self.business_logic.get_collection_version(collection.version_id)
                response["private"][version.collection_id.id] = get_datasets(version)
        return response

    def dataset_migrate(self, collection_id: str, dataset_id: str, dataset_version_id: str) -> Dict[str, str]:
        raw_h5ad_uri = [
            artifact.uri
            for artifact in self.business_logic.get_dataset_artifacts(dataset_version_id)
            if artifact.type == DatasetArtifactType.RAW_H5AD
        ]
        bucket_name, object_key = self.business_logic.s3_provider.parse_s3_uri(raw_h5ad_uri)
        self.business_logic.s3_provider.download_file(bucket_name, object_key, "raw.h5ad")
        cellxgene_schema.migrate("raw.h5ad", "migrated.h5ad", collection_id, dataset_id)
        # TODO: define where this upload bucket actually should be
        upload_bucket_name = os.environ["UPLOAD_BUCKET"]
        dst_uri = f"{dataset_version_id}/migrated.h5ad"
        self.business_logic.s3_provider.upload_file("migrated.h5ad", upload_bucket_name, dst_uri)
        url = f"s3://{upload_bucket_name}/{dst_uri}"
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(collection_id, url, None, dataset_version_id)
        return {"new_dataset_version_id": new_dataset_version_id}

    def collection_migrate(
        self, collection_id: str, datasets: List[Dict[str, str]], can_open_revision: bool
    ) -> List[Dict[str, str]]:
        private_collection_id = collection_id
        if can_open_revision:
            private_collection_id = self.business_logic.create_collection_version(collection_id).version_id.id
        return [
            {
                "collection_id": private_collection_id,
                "dataset_id": dataset["dataset_id"],
                "dataset_version_id": dataset["dataset_version_id"],
            }
            for dataset in datasets
        ]

    def corpus_migrate(self) -> bool:
        self.gather_collections()
        return True

    def migrate(self) -> bool:
        """
        Gets called by the step function at every different step, as defined by `step_name`
        """
        if self.step_name == "corpus_migrate":
            self.corpus_migrate()
        if self.step_name == "collection_migrate":
            collection_id = os.environ["collection_id"]
            datasets = json.loads(os.environ["datasets"])
            can_open_revision = os.environ["can_open_revision"].lower() == "true"
            self.collection_migrate(collection_id, datasets, can_open_revision)
        if self.step_name == "dataset_migrate":
            collection_id = os.environ["collection_id"]
            dataset_id = os.environ["dataset_id"]
            dataset_version_id = os.envion["dataset_version_id"]
            self.dataset_migrate(collection_id, dataset_id, dataset_version_id)
        return True


if __name__ == "__main__":
    business_logic = BusinessLogic(
        DatabaseProvider(),
        None,  # Not required
        None,  # Not required
        S3Provider(),
        UriProvider(),
    )
    step_name = os.environ["STEP_NAME"]
    schema_migrate = SchemaMigrate(business_logic, step_name)
    rv = schema_migrate.migrate()
    sys.exit(rv)
