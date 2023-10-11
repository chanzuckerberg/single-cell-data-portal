import itertools
import json
import logging
import os
import random
from typing import Any, Dict, Iterable, List

from backend.common.corpora_config import CorporaConfig
from backend.common.utils.result_notification import upload_to_slack
from backend.layers.business.business import BusinessLogic
from backend.layers.business.entities import CollectionQueryFilter
from backend.layers.common.entities import (
    CollectionId,
    CollectionVersion,
    CollectionVersionId,
    DatasetArtifactType,
    DatasetProcessingStatus,
    DatasetVersionId,
)
from backend.layers.processing import logger
from backend.layers.processing.process_logic import ProcessingLogic
from backend.layers.thirdparty.schema_validator_provider import SchemaValidatorProvider
from backend.layers.thirdparty.step_function_provider import StepFunctionProvider, sfn_name_generator

logger.configure_logging(level=logging.INFO)


class SchemaMigrate(ProcessingLogic):
    def __init__(self, business_logic: BusinessLogic, schema_validator: SchemaValidatorProvider):
        self.schema_validator = schema_validator
        self.business_logic = business_logic
        self.s3_provider = business_logic.s3_provider  # For compatiblity with ProcessingLogic
        self.artifact_bucket = os.environ.get("ARTIFACT_BUCKET", "test-bucket")
        self.execution_id = os.environ.get("EXECUTION_ID", "test-execution-arn")
        self.logger = logging.getLogger("processing")
        self.local_path: str = "."  # Used for testing
        self.limit_migration = os.environ.get("LIMIT_MIGRATION", False)  # Run a small migration for testing
        self.limit_select = 2  # Number of collections to migrate
        self.schema_version = schema_validator.get_current_schema_version()

    def limit_collections(self) -> Iterable[CollectionVersion]:
        published_collections = [*self.business_logic.get_collections(CollectionQueryFilter(is_published=True))]
        unpublished_collections = [*self.business_logic.get_collections(CollectionQueryFilter(is_published=False))]
        if self.limit_migration:
            select = self.limit_select // 2
            if len(unpublished_collections) >= select:
                unpublished_collections = random.sample(unpublished_collections, self.limit_select // 2)
            if len(published_collections) >= select:
                published_collections = random.sample(published_collections, self.limit_select // 2)
        return itertools.chain(unpublished_collections, published_collections)

    def gather_collections(self, auto_publish: bool) -> List[Dict[str, str]]:
        """
        This function is used to gather all the collections and their datasets that will be migrated
        :return: A dictionary with the following structure:
        [
            {"can_publish": "true", "collection_id": "<collection_id>", "collection_version_id":
            "<collection_version_id>"},
            {"can_publish": "false", "collection_id": "<collection_id>", "collection_version_id":
            "<collection_version_id>"}
            ...
        ]

        :param auto_publish: bool - if False, coerce can_publish to False for all collections. if True, determine
        can_publish on collection-by-collection basis based on business logic
        """
        response = []

        has_revision = set()
        # iterates over unpublished collections first, so published versions are skipped if there is an active revision
        for collection in self.limit_collections():
            _resp = {}
            if collection.is_published() and collection.collection_id.id in has_revision:
                continue

            if not auto_publish:
                # auto_publish is off for this migration
                _resp["can_publish"] = str(False)
            elif collection.is_published():
                # published collection without an active revision
                _resp["can_publish"] = str(True)
            elif collection.is_unpublished_version():
                # active revision of a published collection.
                has_revision.add(collection.collection_id.id)  # revision found, skip published version
                _resp["can_publish"] = str(False)
            elif collection.is_initial_unpublished_version():
                # unpublished collection
                _resp["can_publish"] = str(False)
            _resp.update(
                collection_id=collection.collection_id.id,
                collection_version_id=collection.version_id.id,
            )
            response.append(_resp)
        return response

    def dataset_migrate(
        self, collection_version_id: str, collection_id: str, dataset_id: str, dataset_version_id: str
    ) -> Dict[str, str]:
        raw_h5ad_uri = [
            artifact.uri
            for artifact in self.business_logic.get_dataset_artifacts(DatasetVersionId(dataset_version_id))
            if artifact.type == DatasetArtifactType.RAW_H5AD
        ][0]
        source_bucket_name, source_object_key = self.s3_provider.parse_s3_uri(raw_h5ad_uri)
        self.s3_provider.download_file(source_bucket_name, source_object_key, "previous_schema.h5ad")
        migrated_file = "migrated.h5ad"
        self.schema_validator.migrate("previous_schema.h5ad", migrated_file, collection_id, dataset_id)
        key_prefix = self.get_key_prefix(dataset_version_id)
        uri = self.upload_artifact(migrated_file, key_prefix, self.artifact_bucket)
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            CollectionVersionId(collection_version_id),
            uri,
            file_size=0,  # TODO: this shouldn't be needed but it gets around a 404 for HeadObject
            existing_dataset_version_id=DatasetVersionId(dataset_version_id),
            start_step_function=False,  # The schema_migration sfn will start the ingest sfn
        )
        sfn_name = sfn_name_generator(dataset_version_id, prefix="migrate")
        return {
            "collection_version_id": collection_version_id,
            "dataset_version_id": new_dataset_version_id.id,
            "uri": uri,
            "sfn_name": sfn_name,
        }

    def collection_migrate(self, collection_id: str, collection_version_id: str, can_publish: bool) -> Dict[str, Any]:
        # Get datasets from collection
        version = self.business_logic.get_collection_version(CollectionVersionId(collection_version_id))
        datasets = [dataset for dataset in version.datasets if not self.check_dataset_is_latest_schema_version(dataset)]
        # Filter out datasets that are already on the current schema version
        if not datasets:
            # Handles the case were the collection has no datasets or all datasets are already migrated.
            if len(version.datasets) == 0:
                self.logger.info("Collection has no datasets")
            else:
                self.logger.info(
                    "All datasets in the collection have been migrated", extra={"dataset_count": len(version.datasets)}
                )
            response = {
                "can_publish": str(False),  # skip publishing, because the collection is already published and no
                # revision is created, or the collection is private or a revision.
                "collection_version_id": collection_version_id,
                "datasets": [],
                "no_datasets": str(True),
            }
        else:
            if version.is_published():
                # Create a new collection version(revision) if the collection is already published
                private_collection_version_id = self.business_logic.create_collection_version(
                    CollectionId(collection_id)
                ).version_id.id
            else:
                private_collection_version_id = collection_version_id

            response = {
                "can_publish": str(can_publish),
                "collection_version_id": private_collection_version_id,
                # ^^^ The top level fields are used for handling error cases in the AWS SFN.
                "datasets": [
                    {
                        "can_publish": str(can_publish),
                        "collection_id": collection_id,
                        "collection_version_id": private_collection_version_id,
                        "dataset_id": dataset.dataset_id.id,
                        "dataset_version_id": dataset.version_id.id,
                    }
                    for dataset in datasets
                    if dataset.status.processing_status == DatasetProcessingStatus.SUCCESS
                    # Filter out datasets that are not successfully processed
                ]
                # The repeated fields in datasets is required for the AWS SFN job that uses it.
            }

            if not response["datasets"]:
                # Handles the case were the collection has no processed datasets
                response["no_datasets"] = str(True)
        self._store_sfn_response("publish_and_cleanup", version.collection_id.id, response)
        return response

    def publish_and_cleanup(self, collection_version_id: str, can_publish: bool) -> list:
        errors = []
        collection_version = self.business_logic.get_collection_version(CollectionVersionId(collection_version_id))
        object_keys_to_delete = []

        # Get the datasets that were processed
        extra_info = self._retrieve_sfn_response("publish_and_cleanup", collection_version.collection_id.id)
        processed_datasets = {d["dataset_id"]: d["dataset_version_id"] for d in extra_info["datasets"]}

        # Process datasets errors
        for dataset in collection_version.datasets:
            dataset_id = dataset.dataset_id.id
            dataset_version_id = dataset.version_id.id
            _log_extras = {
                "dataset_id": dataset_id,
                "dataset_version_id": dataset_version_id,
            }
            if dataset_id not in processed_datasets:
                self.logger.info("skipping dataset", extra=_log_extras)
                continue
            # filepath to clean-up uses dataset_version_id from the replaced version; accessing with dataset_id as key
            previous_dataset_version_id = processed_datasets[dataset_id]
            _log_extras["previous_dataset_version_id"] = previous_dataset_version_id
            self.logger.info("checking dataset", extra=_log_extras)
            key_prefix = self.get_key_prefix(previous_dataset_version_id)
            object_keys_to_delete.append(f"{key_prefix}/migrated.h5ad")
            if not self.check_dataset_is_latest_schema_version(dataset):
                error = {
                    "message": "Did Not Migrate.",
                    "collection_id": collection_version.collection_id.id,
                    "collection_version_id": collection_version_id,
                    "dataset_version_id": dataset_version_id,
                    "dataset_id": dataset_id,
                    "rollback": False,
                }
                self.logger.error(error)
                errors.append(error)
            elif dataset.status.processing_status != DatasetProcessingStatus.SUCCESS:
                error = {
                    "message": dataset.status.validation_message,
                    "dataset_status": dataset.status.to_dict(),
                    "collection_id": collection_version.collection_id.id,
                    "collection_version_id": collection_version_id,
                    "dataset_version_id": dataset_version_id,
                    "dataset_id": dataset_id,
                    "rollback": True,
                }
                self.logger.error(error)
                errors.append(error)
            else:
                self.logger.info("checked dataset")
        self.logger.info(
            "Deleting files", extra={"artifact_bucket": self.artifact_bucket, "object_keys": object_keys_to_delete}
        )
        self.s3_provider.delete_files(self.artifact_bucket, object_keys_to_delete)
        if errors:
            self._store_sfn_response("report", collection_version_id, errors)
        elif can_publish:
            self.business_logic.publish_collection_version(collection_version.version_id)
        return errors

    def _store_sfn_response(self, step_name, file_name, response: Dict[str, str]):
        """

        :param step_name: The step that will use this file
        :param file_name: a unique name to describe this job
        :param response: the response to store as json.
        """
        file_name = f"{file_name}.json"
        local_file = os.path.join(self.local_path, file_name)
        key_name = self.get_key_prefix(f"schema_migration/{self.execution_id}/{step_name}/{file_name}")
        with open(local_file, "w") as f:
            json.dump(response, f)
        self.s3_provider.upload_file(local_file, self.artifact_bucket, key_name, {})
        self.logger.info(
            "Uploaded to S3", extra={"file_name": local_file, "bucket": self.artifact_bucket, "key": key_name}
        )

    def _retrieve_sfn_response(self, step_name, file_name):
        """
        retrieve the JSON responses to be used by this step
        :param step_name: the step that the response is intended for
        :param file_name: a unique name of the file.
        :return: the contents of the JSON file
        """
        file_name = f"{file_name}.json"
        local_file = os.path.join(self.local_path, "data.json")
        key_name = self.get_key_prefix(f"schema_migration/{self.execution_id}/{step_name}/{file_name}")
        self.s3_provider.download_file(self.artifact_bucket, key_name, local_file)
        with open(local_file, "r") as f:
            data = json.load(f)
        self.s3_provider.delete_files(self.artifact_bucket, [key_name])  # delete after reading.
        return data

    def error_wrapper(self, func, file_name: str):
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                self.logger.exception(f"Error in {func.__name__}", extra={"input": {"args": args, "kwargs": kwargs}})
                self._store_sfn_response(
                    "report", file_name, {"step": func.__name__, "error": str(e), "args": args, "kwargs": kwargs}
                )
                raise e

        return wrapper

    def report(self) -> str:
        try:
            report = dict(errors=[])
            error_s3_keys = list(
                self.s3_provider.list_directory(
                    self.artifact_bucket, self.get_key_prefix(f"schema_migration/{self.execution_id}/report")
                )
            )
            self.logger.info("Error files found", extra={"error_files": len(error_s3_keys)})
            for s3_key in error_s3_keys:
                local_file = os.path.join(self.local_path, "data.json")
                self.s3_provider.download_file(self.artifact_bucket, s3_key, local_file)
                with open(local_file, "r") as f:
                    jn = json.load(f)
                report["errors"].append(jn)
            self.logger.info("Report", extra=report)
            report_str = json.dumps(report, indent=4, sort_keys=True)

            # Cleanup S3 files
            self.s3_provider.delete_files(self.artifact_bucket, error_s3_keys)
            report_message = "Schema migration results."
            # if report["errors"]:
            #     report_message += " @sc-oncall-eng"
            self._upload_to_slack("schema_migration_report.json", report_str, report_message)
            return report
        except Exception as e:
            self.logger.exception("Failed to generate report")
            raise e

    def _upload_to_slack(self, filename: str, contents, initial_comment: str) -> None:
        slack_token = CorporaConfig().slack_reporter_secret
        slack_channel = CorporaConfig().slack_reporter_channel
        response = upload_to_slack(filename, contents, initial_comment, slack_channel, slack_token)
        self.logger.info("Uploaded to slack", extra={"response": response})

    def migrate(self, step_name) -> bool:
        """
        Gets called by the step function at every different step, as defined by `step_name`
        """
        self.logger.info(f"Starting {step_name}", extra={"step": step_name})
        if step_name == "gather_collections":
            gather_collections = self.error_wrapper(self.gather_collections, "gather_collections")
            auto_publish = os.environ["AUTO_PUBLISH"].lower() == "true"
            response = gather_collections(auto_publish)
        elif step_name == "collection_migrate":
            collection_id = os.environ["COLLECTION_ID"]
            collection_version_id = os.environ["COLLECTION_VERSION_ID"]
            can_publish = os.environ["CAN_PUBLISH"].lower() == "true"
            collection_migrate = self.error_wrapper(self.collection_migrate, collection_id)
            response = collection_migrate(
                collection_id=collection_id,
                collection_version_id=collection_version_id,
                can_publish=can_publish,
            )
        elif step_name == "dataset_migrate":
            collection_version_id = os.environ["COLLECTION_VERSION_ID"]
            collection_id = os.environ["COLLECTION_ID"]
            dataset_id = os.environ["DATASET_ID"]
            dataset_version_id = os.environ["DATASET_VERSION_ID"]
            dataset_migrate = self.error_wrapper(self.dataset_migrate, f"{collection_version_id}_{dataset_id}")
            response = dataset_migrate(
                collection_version_id=collection_version_id,
                collection_id=collection_id,
                dataset_id=dataset_id,
                dataset_version_id=dataset_version_id,
            )
        elif step_name == "collection_publish":
            collection_version_id = os.environ["COLLECTION_VERSION_ID"]
            can_publish = os.environ["CAN_PUBLISH"].lower() == "true"
            publish_and_cleanup = self.error_wrapper(self.publish_and_cleanup, collection_version_id)
            response = publish_and_cleanup(collection_version_id=collection_version_id, can_publish=can_publish)
        elif step_name == "report":
            response = self.report()
        self.logger.info("output", extra={"response": response})
        sfn_client = StepFunctionProvider().client
        sfn_client.send_task_success(taskToken=os.environ["TASK_TOKEN"], output=json.dumps(response))
        return True
