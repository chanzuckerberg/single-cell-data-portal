import itertools
import json
import logging
import os
from typing import Dict, List

import cellxgene_schema

from backend.layers.business.business import BusinessLogic
from backend.layers.business.entities import CollectionQueryFilter
from backend.layers.common.entities import (
    CollectionId,
    CollectionVersionId,
    DatasetArtifactType,
    DatasetProcessingStatus,
    DatasetVersionId,
)
from backend.layers.processing import logger

logger.configure_logging(level=logging.INFO)


class SchemaMigrate:
    def __init__(self, business_logic: BusinessLogic):
        self.business_logic = business_logic

    def gather_collections(self) -> Dict[str, Dict[str, List[Dict[str, str]]]]:
        """
        This function is used to gather all the collections and their datasets that will be migrated
        :return: A dictionary with the following structure:
        {
            "published": {
                "<collection_id>":
                    [
                        {"dataset_id": "<dataset_id>", "dataset_version_id": "<dataset_version_id>"},
                        {...},
                         ...
                     ],
                ...
            },
            "revision": {
                "<collection_version_id>":
                    [
                        {"dataset_id": "<dataset_id>", "dataset_version_id": "<dataset_version_id>"},
                        {...},
                         ...
                     ],
                ...
            },
            "private": {
                "<collection_id>":
                    [
                        {"dataset_id": "<dataset_id>", "dataset_version_id": "<dataset_version_id>"},
                        {...},
                         ...
                     ],
                ...
            },
        }
        """

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
            for artifact in self.business_logic.get_dataset_artifacts(DatasetVersionId(dataset_version_id))
            if artifact.type == DatasetArtifactType.RAW_H5AD
        ]
        bucket_name, object_key = self.business_logic.s3_provider.parse_s3_uri(raw_h5ad_uri)
        self.business_logic.s3_provider.download_file(bucket_name, object_key, "previous_schema.h5ad")
        cellxgene_schema.migrate("previous_schema.h5ad", "migrated.h5ad", collection_id, dataset_id)
        upload_bucket = os.environ["ARTIFACT_BUCKET"]
        dst_uri = f"{dataset_version_id}/migrated.h5ad"
        self.business_logic.s3_provider.upload_file("migrated.h5ad", upload_bucket, dst_uri)
        url = f"s3://{upload_bucket}/{dst_uri}"
        return {"collection_id": collection_id, "dataset_version_id": dataset_version_id, "url": url}

    def collection_migrate(
        self, collection_id: str, datasets: List[Dict[str, str]], can_open_revision: bool
    ) -> List[Dict[str, str]]:
        private_collection_id = collection_id
        if can_open_revision:
            private_collection_id = self.business_logic.create_collection_version(
                CollectionId(collection_id)
            ).version_id.id
        return [
            {
                "collection_id": private_collection_id,
                "dataset_id": dataset["dataset_id"],
                "dataset_version_id": dataset["dataset_version_id"],
            }
            for dataset in datasets
        ]

    def publish_and_cleanup(self, collection_id: str, can_publish: bool) -> Dict[str, str]:
        errors = dict()
        collection_version_id = CollectionVersionId(collection_id)
        collection_version = self.business_logic.get_collection_version(collection_version_id)
        current_schema_version = cellxgene_schema.schema.get_current_schema_version()
        object_keys_to_delete = []
        for dataset in collection_version.datasets:
            dataset_version_id = dataset.version_id.id
            object_keys_to_delete.append(f"{dataset_version_id}/migrated.h5ad")
            if dataset.metadata.schema_version != current_schema_version:
                errors[dataset_version_id] = "Did Not Migrate."
            elif dataset.status.processing_status != DatasetProcessingStatus.SUCCESS:
                errors[dataset_version_id] = dataset.status.validation_message

        artifact_bucket = os.environ["ARTIFACT_BUCKET"]
        if errors:
            self._store_in_s3("report", collection_id, errors)
        elif can_publish:
            self.business_logic.publish_collection_version(collection_version_id)
        self.business_logic.s3_provider.delete_files(artifact_bucket, object_keys_to_delete)
        return errors

    def _store_in_s3(self, step_name, file_name, response: Dict[str, str]):
        """

        :param step_name: The step that will use this file
        :param file_name: a unique name to describe this job
        :param response: the response to store as json.
        """
        execution_arn = os.environ["EXECUTION_ARN"]
        with open("response.json", "w") as f:
            json.dump(response, f)
        self.business_logic.s3_provider.upload_file(
            "response.json",
            os.environ["ARTIFACT_BUCKET"],
            f"schema_migration/{execution_arn}/{step_name}/{file_name}.json",
        )

    def report(self):
        bucket = os.environ["ARTIFACT_BUCKET"]
        error_files = list(
            self.business_logic.s3_provider.list_directory(
                bucket, f"schema_migration/" f"{os.environ['EXECUTION_ARN']}/report"
            )
        )
        for file in error_files:
            self.business_logic.s3_provider.download_file(bucket, file, file)
            with open(file, "r") as f:
                json.load(f)
            pass  # TODO generate a report from the errors.
        self.business_logic.s3_provider.delete_files(bucket, error_files)

    def migrate(self, step_name) -> bool:
        """
        Gets called by the step function at every different step, as defined by `step_name`
        """
        if step_name == "gather_collections":
            self.gather_collections()
        if step_name == "collection_migrate":
            collection_id = os.environ["collection_id"]
            datasets = json.loads(os.environ["datasets"])
            can_open_revision = os.environ["can_open_revision"].lower() == "true"
            self.collection_migrate(collection_id, datasets, can_open_revision)
        if step_name == "dataset_migrate":
            collection_id = os.environ["collection_id"]
            dataset_id = os.environ["dataset_id"]
            dataset_version_id = os.environ["dataset_version_id"]
            self.dataset_migrate(collection_id, dataset_id, dataset_version_id)
        if step_name == "publish_and_cleanup":
            collection_id = os.environ["collection_id"]
            can_publish = os.environ["can_publish"].lower() == "true"
            self.publish_and_cleanup(collection_id, can_publish)
        if step_name == "report":
            self.report()
        return True
