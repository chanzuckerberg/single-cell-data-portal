import itertools
import json
import logging
import os
from typing import Dict, List

from cellxgene_schema import schema

from backend.common.corpora_config import CorporaConfig
from backend.common.utils.result_notification import upload_to_slack
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
from backend.layers.thirdparty.step_function_provider import StepFunctionProvider

logger.configure_logging(level=logging.INFO)


class SchemaMigrate:
    def __init__(self, business_logic: BusinessLogic):
        self.business_logic = business_logic
        self.bucket = os.environ.get("ARTIFACT_BUCKET", "test-bucket")
        self.execution_arn = os.environ.get("EXECUTION_ARN", "test-execution-arn")
        self.logger = logging.getLogger(__name__)

    def gather_collections(self) -> List[Dict[str, str]]:
        """
        This function is used to gather all the collections and their datasets that will be migrated
        :return: A dictionary with the following structure:
        [
            {"can_open_revision": True, "collection_id": "<collection_id>"},
            {"can_open_revision": False, "collection_id": "<collection_id>", "collection_version_id":
            "<collection_version_id>"},
            {"can_open_revision": False, "collection_id": "<collection_id>"},
            ...
        ]
        """

        published_collections = self.business_logic.get_collections(CollectionQueryFilter(is_published=True))
        unpublished_collections = self.business_logic.get_collections(CollectionQueryFilter(is_published=False))

        response = []

        collections = itertools.chain(unpublished_collections, published_collections)
        # evaluate unpublished collections first, so that published versions are skipped if there is an active revision
        has_revision = []  # list of collections to skip if published with an active revision
        for collection in collections:

            # TODO <testing code>
            if not collection.metadata.name.startswith("TestSchemaMigrate"):
                print("Skipping collection")
                continue
            # TODO </testing code>

            if collection.is_published() and collection.collection_id not in has_revision:
                # published collection without an active revision
                response.append(
                    dict(
                        can_open_revision="True",
                        collection_id=collection.collection_id.id,
                        collection_version_id=collection.version_id.id,
                    )
                )
            elif collection.is_unpublished_version():
                # published collection with an active revision
                has_revision.append(collection.collection_id)  # revision found, skip published version
                response.append(
                    dict(
                        can_open_revision="False",
                        collection_id=collection.collection_id.id,
                        collection_version_id=collection.version_id.id,
                    )
                )
                # include collection version id
            elif collection.is_initial_unpublished_version():
                # unpublished collection
                response.append(
                    dict(
                        can_open_revision="False",
                        collection_id=collection.collection_id.id,
                        collection_version_id=collection.version_id.id,
                    )
                )
        return response

    def dataset_migrate(self, collection_id: str, dataset_id: str, dataset_version_id: str) -> Dict[str, str]:
        raw_h5ad_uri = [
            artifact.uri
            for artifact in self.business_logic.get_dataset_artifacts(DatasetVersionId(dataset_version_id))
            if artifact.type == DatasetArtifactType.RAW_H5AD
        ][0]
        bucket_name, object_key = self.business_logic.s3_provider.parse_s3_uri(raw_h5ad_uri)
        self.business_logic.s3_provider.download_file(bucket_name, object_key, "previous_schema.h5ad")
        schema.migrate("previous_schema.h5ad", "migrated.h5ad", collection_id, dataset_id)
        upload_bucket = os.environ["ARTIFACT_BUCKET"]
        dst_uri = f"{dataset_version_id}/migrated.h5ad"
        self.business_logic.s3_provider.upload_file("migrated.h5ad", upload_bucket, dst_uri, {})
        url = f"s3://{upload_bucket}/{dst_uri}"
        return {"collection_id": collection_id, "dataset_version_id": dataset_version_id, "url": url}

    def collection_migrate(
        self, collection_id: str, collection_version_id: str, can_open_revision: bool
    ) -> List[Dict[str, str]]:
        private_collection_id = collection_version_id

        # Get datasets from collection
        version = self.business_logic.get_collection_version(CollectionVersionId(collection_version_id))
        datasets = []
        for dataset in version.datasets:
            datasets.append({"dataset_id": dataset.dataset_id.id, "dataset_version_id": dataset.version_id.id})

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
        current_schema_version = schema.get_current_schema_version()
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
        file_name = f"{file_name}.json"
        with open(file_name, "w") as f:
            json.dump(response, f)
        self.business_logic.s3_provider.upload_file(
            file_name,
            self.bucket,
            f"schema_migration/{self.execution_arn}/{step_name}/{file_name}",
        )

    def error_decorator(self, func, file_name: str):
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                self.logger.exception(f"Error in {func.__name__}", extra={"input": {"args": args, "kwargs": kwargs}})
                self._store_in_s3(
                    func.__name__, file_name, {"step": func.__name__, "error": str(e), "args": args, "kwargs": kwargs}
                )
                raise e

        return wrapper

    def report(self, local_path: str = ".") -> str:
        try:
            report = dict(errors=[])
            error_files = list(
                self.business_logic.s3_provider.list_directory(
                    self.bucket, f"schema_migration/{self.execution_arn}/report"
                )
            )
            for file in error_files:
                local_file = os.path.join(local_path, file)
                self.business_logic.s3_provider.download_file(self.bucket, file, local_file)
                with open(local_file, "r") as f:
                    jn = json.load(f)
                report["errors"].append(jn)
            self.logger.info("Report", extra=report)
            report_str = json.dumps(report)
            self.business_logic.s3_provider.delete_files(self.bucket, error_files)
            self._upload_to_slack("schema_migration_report.json", report_str, "Schema migration results.")
            return report
        except Exception as e:
            self.logger.exception("Failed to generate report")
            raise e

    def _upload_to_slack(self, filename: str, contents, initial_comment: str) -> None:
        slack_token = CorporaConfig().slack_reporter_secret
        slack_channel = CorporaConfig().slack_reporter_channel
        response = upload_to_slack(filename, contents, initial_comment, slack_channel, slack_token)
        return response

    def migrate(self, step_name) -> bool:
        """
        Gets called by the step function at every different step, as defined by `step_name`
        """
        self.logger.info(f"Starting {step_name}", extra={"step": step_name})
        if step_name == "gather_collections":
            gather_collections = self.error_decorator(self.gather_collections, "gather_collections")
            response = gather_collections()
        elif step_name == "collection_migrate":
            collection_id = os.environ["COLLECTION_ID"]
            collection_version_id = os.environ["COLLECTION_VERSION_ID"]
            can_open_revision = os.environ["CAN_OPEN_REVISION"].lower() == "true"
            collection_migrate = self.error_decorator(self.collection_migrate, collection_id)
            response = collection_migrate(
                collection_id=collection_id,
                collection_version_id=collection_version_id,
                can_open_revision=can_open_revision,
            )
        elif step_name == "dataset_migrate":
            collection_id = os.environ["COLLECTION_ID"]
            dataset_id = os.environ["DATASET_ID"]
            dataset_version_id = os.environ["DATASET_VERSION_ID"]
            dataset_migrate = self.error_decorator(self.dataset_migrate, f"{collection_id}_{dataset_id}")
            response = dataset_migrate(
                collection_id=collection_id, dataset_id=dataset_id, dataset_version_id=dataset_version_id
            )
        elif step_name == "collection_publish":
            collection_id = os.environ["COLLECTION_ID"]
            can_publish = os.environ["CAN_PUBLISH"].lower() == "true"
            publish_and_cleanup = self.error_decorator(self.publish_and_cleanup, collection_id)
            response = publish_and_cleanup(collection_id=collection_id, can_publish=can_publish)
        elif step_name == "report":
            response = self.report()
        self.logger.info("Response", extra=response)
        sfn_client = StepFunctionProvider().client
        sfn_client.send_task_success(taskToken=os.environ["TASK_TOKEN"], output=json.dumps(response))
        return True
