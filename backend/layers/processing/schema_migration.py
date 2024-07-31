import itertools
import json
import logging
import os
import random
from typing import Dict, Iterable, List, Tuple

from backend.common.corpora_config import CorporaConfig
from backend.common.utils.json import CustomJSONEncoder
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
        self.artifact_bucket = os.environ.get("ARTIFACT_BUCKET", "artifact-bucket")
        self.execution_id = os.environ.get("EXECUTION_ID", "test-execution-arn")
        self.logger = logging.getLogger("processing")
        self.local_path: str = "."  # Used for testing
        self.limit_migration = os.environ.get("LIMIT_MIGRATION", 0)  # Run a small migration for testing
        self._schema_version = None

    @property
    def schema_version(self):
        if not self._schema_version:
            self._schema_version = self.schema_validator.get_current_schema_version()
        return self._schema_version

    def fetch_collections(self) -> Iterable[CollectionVersion]:
        published_collections = [*self.business_logic.get_collections(CollectionQueryFilter(is_published=True))]
        unpublished_collections = [*self.business_logic.get_collections(CollectionQueryFilter(is_published=False))]
        return itertools.chain(unpublished_collections, published_collections)

    def gather_collections(self) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        This function is used to gather all the collections and their datasets that will be migrated
        A json file is created and uploaded to S3 with the list of collections and datasets that will be migrated. It
        has the following structure:
        [
            {"collection_id": "<collection_id>", "collection_version_id": "<collection_version_id>"},
            {"collection_id": "<collection_id>", "collection_version_id": "<collection_version_id>"}
            ...
        ]

        :return: the response returned to the step function and the list of collections to be migrated
        """
        response_for_span_collections = []

        has_migration_revision = set()
        # iterates over unpublished collections first, so published versions are skipped if there is an active revision
        for collection in self.fetch_collections():
            if collection.is_published() and collection.collection_id.id in has_migration_revision:
                continue

            if collection.is_auto_version:
                has_migration_revision.add(collection.collection_id.id)  # migration revision found, skip published

            _resp = {
                "collection_id": collection.collection_id.id,
                "collection_version_id": collection.version_id.id,
                "execution_id": self.execution_id,
            }
            response_for_span_collections.append(_resp)

        # For testing purposes, only migrate a randomly sampled subset of the collections gathered
        limit = int(self.limit_migration) if isinstance(self.limit_migration, str) else self.limit_migration
        if limit > 0:
            response_for_span_collections = random.sample(response_for_span_collections, limit)
        key_name = self._store_sfn_response("span_collections", "collections", response_for_span_collections)
        response_for_sfn = {"key_name": key_name}
        return response_for_sfn, response_for_span_collections

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
        reported_changes = self.schema_validator.migrate(
            "previous_schema.h5ad", migrated_file, collection_id, dataset_id
        )
        if reported_changes:
            self._store_sfn_response(
                "report/migrate_changes",
                f"{collection_id}_{dataset_id}",
                {f"{collection_id}_{dataset_id}": reported_changes},
            )
        key_prefix = self.get_key_prefix(dataset_version_id)
        uri = self.upload_artifact(migrated_file, key_prefix, self.artifact_bucket)
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            CollectionVersionId(collection_version_id),
            uri,
            file_size=0,  # TODO: this shouldn't be needed but it gets around a 404 for HeadObject
            current_dataset_version_id=DatasetVersionId(dataset_version_id),
            start_step_function=False,  # The schema_migration sfn will start the ingest sfn
        )
        sfn_name = sfn_name_generator(new_dataset_version_id, prefix="migrate")
        return {
            "collection_version_id": collection_version_id,
            "dataset_version_id": new_dataset_version_id.id,
            "uri": uri,
            "sfn_name": sfn_name,
            "execution_id": self.execution_id,
        }

    def collection_migrate(
        self,
        collection_id: str,
        collection_version_id: str,
    ) -> Tuple[Dict[str, str], Dict[str, str], List[Dict[str, str]]]:
        """
        This function is used to migrate a collection and its datasets to the latest schema version.

        :param collection_id: the canonical collection id
        :param collection_version_id: the collection version to migrate
        :return: the response retuned to the step function, the response for the log_errors_and_cleanup step function,
        and the list of datasets to be migrated
        """
        # Get datasets from collection
        version = self.business_logic.get_collection_version(CollectionVersionId(collection_version_id))
        # Filter out datasets that are already on the current schema version
        datasets = [dataset for dataset in version.datasets if not self.check_dataset_is_latest_schema_version(dataset)]

        # Generate canonical collection url
        collection_url = self.business_logic.get_collection_url(version.collection_id.id)

        if not datasets:
            # Handles the case were the collection has no datasets or all datasets are already migrated.
            if len(version.datasets) == 0:
                self.logger.info("Collection has no datasets")
            else:
                self.logger.info(
                    "All datasets in the collection have been migrated", extra={"dataset_count": len(version.datasets)}
                )
            response_for_dataset_migrate = []
            response_for_log_errors_and_cleanup = {
                "collection_version_id": collection_version_id,
            }
            response_for_sfn = {"collection_version_id": collection_version_id}
        else:
            if version.is_published():
                # Create a new collection version(revision) if the collection is already published
                private_collection_version_id = self.business_logic.create_collection_version(
                    CollectionId(collection_id),
                    is_auto_version=True,
                ).version_id.id
            else:
                private_collection_version_id = collection_version_id
            response_for_dataset_migrate = [
                {
                    "collection_id": collection_id,
                    "collection_url": collection_url,
                    "collection_version_id": private_collection_version_id,
                    "dataset_id": dataset.dataset_id.id,
                    "dataset_version_id": dataset.version_id.id,
                    "execution_id": self.execution_id,
                }
                for dataset in datasets
                if dataset.status.processing_status == DatasetProcessingStatus.SUCCESS
                # Filter out datasets that are not successfully processed
            ]
            response_for_log_errors_and_cleanup = {
                "collection_version_id": private_collection_version_id,
            }
            response_for_sfn = {"collection_version_id": private_collection_version_id}
        response_for_log_errors_and_cleanup["datasets"] = response_for_dataset_migrate
        response_for_log_errors_and_cleanup["collection_url"] = collection_url

        response_for_sfn["execution_id"] = self.execution_id

        self._store_sfn_response(
            "log_errors_and_cleanup", version.collection_id.id, response_for_log_errors_and_cleanup
        )

        if response_for_dataset_migrate:
            key_name = self._store_sfn_response("span_datasets", version.collection_id.id, response_for_dataset_migrate)
            response_for_sfn["key_name"] = key_name
        return (response_for_sfn, response_for_log_errors_and_cleanup, response_for_dataset_migrate)

    def log_errors_and_cleanup(self, collection_version_id: str) -> list:
        errors = []
        rolled_back_datasets = []
        collection_version = self.business_logic.get_collection_version(CollectionVersionId(collection_version_id))
        object_keys_to_delete = []

        # Get the datasets that were processed
        extra_info = self._retrieve_sfn_response("log_errors_and_cleanup", collection_version.collection_id.id)
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
            if dataset.status.processing_status != DatasetProcessingStatus.SUCCESS:
                error = {
                    "message": dataset.status.validation_message,
                    "dataset_status": dataset.status.to_dict(),
                    "collection_id": collection_version.collection_id.id,
                    "collection_version_id": collection_version_id,
                    "dataset_version_id": dataset_version_id,
                    "dataset_id": dataset_id,
                    "rollback": True,
                }
                # If the dataset is not successfully processed, rollback to the version from before migration
                self.business_logic.restore_previous_dataset_version(
                    CollectionVersionId(collection_version_id), dataset.dataset_id
                )
                rolled_back_datasets.append(dataset)
                self.logger.error(error)
                errors.append(error)
            elif not self.check_dataset_is_latest_schema_version(dataset):
                error_message = "Did Not Migrate"
                dataset_status = "n/a"
                if dataset.status is not None:
                    error_message = dataset.status.validation_message
                    dataset_status = dataset.status.to_dict()
                error = {
                    "message": error_message,
                    "dataset_status": dataset_status,
                    "collection_id": collection_version.collection_id.id,
                    "collection_version_id": collection_version_id,
                    "dataset_version_id": dataset_version_id,
                    "dataset_id": dataset_id,
                    "rollback": False,
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
            self._store_sfn_response("report/errors", collection_version_id, errors)
            # clean up artifacts for any now-orphaned, rolled back datasets
            if rolled_back_datasets:
                # TODO: replace with async job to delete orphaned dataset version DB rows + artifacts
                self.business_logic.delete_dataset_versions(rolled_back_datasets)
        return errors

    def _store_sfn_response(self, directory: str, file_name: str, response: Dict[str, str]):
        """
        :param directory: The subdirectory in which to store the response,
        :param file_name: a unique name to describe this job
        :param response: the response to store as json.
        """
        file_name = f"{file_name}.json"
        local_file = os.path.join(self.local_path, file_name)
        key_name = self.get_key_prefix(f"schema_migration/{self.execution_id}/{directory}/{file_name}")
        with open(local_file, "w") as f:
            json.dump(response, f, cls=CustomJSONEncoder)
        self.s3_provider.upload_file(local_file, self.artifact_bucket, key_name, {})
        self.logger.info(
            "Uploaded to S3", extra={"file_name": local_file, "bucket": self.artifact_bucket, "key": key_name}
        )
        return key_name

    def _retrieve_sfn_response(self, directory: str, file_name: str):
        """
        retrieve the JSON responses to be used by this step
        :param directory: the subdirectory from which to retrieve the response
        :param file_name: the filename to retrieve.
        :return: the contents of the JSON file
        """
        file_name = f"{file_name}.json"
        local_file = os.path.join(self.local_path, "data.json")
        key_name = self.get_key_prefix(f"schema_migration/{self.execution_id}/{directory}/{file_name}")
        self.s3_provider.download_file(self.artifact_bucket, key_name, local_file)
        with open(local_file, "r") as f:
            data = json.load(f)
        self.logger.info(
            "Downloaded from S3",
            extra={"file_name": local_file, "bucket": self.artifact_bucket, "key": key_name, "data": data},
        )
        self.s3_provider.delete_files(self.artifact_bucket, [key_name])  # delete after reading.
        return data

    def error_wrapper(self, func, file_name: str):
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                self.logger.exception(f"Error in {func.__name__}", extra={"input": {"args": args, "kwargs": kwargs}})
                self._store_sfn_response(
                    "report/errors", file_name, {"step": func.__name__, "error": str(e), "args": args, "kwargs": kwargs}
                )
                raise e

        return wrapper

    def report(self, artifact_bucket=None, execution_id=None, dry_run=False) -> dict:
        """
        Generate a report of the schema migration process. This function will download all the error and migration
        :param artifact_bucket: The bucket where the schema migration artifacts are stored.
        :param execution_id: the execution id of the AWS SFN schema migration in progress.
        :param dry_run: If dry_run is True, then a report will be returned without deleting any s3 assets or report to
            slack.
        :return: a json report of the schema migration process
        """
        artifact_bucket = artifact_bucket or self.artifact_bucket
        execution_id = execution_id or self.execution_id

        try:
            report = dict(errors=[], migrate_changes=[])

            def retrieve_report_files_from_s3(message_type: str):
                s3_keys = list(
                    self.s3_provider.list_directory(
                        artifact_bucket,
                        self.get_key_prefix(f"schema_migration/{execution_id}/report/{message_type}"),
                    )
                )
                self.logger.info("Subdirectory Count", extra={"message_type": message_type, "count": len(s3_keys)})
                for s3_key in s3_keys:
                    local_file = os.path.join(self.local_path, "data.json")
                    self.s3_provider.download_file(artifact_bucket, s3_key, local_file)
                    with open(local_file, "r") as f:
                        jn = json.load(f)
                    report[message_type].append(jn)

            retrieve_report_files_from_s3("errors")
            retrieve_report_files_from_s3("migrate_changes")
            if not dry_run:
                self.logger.info("Report", extra=report)
                report_str = json.dumps(report, indent=4, sort_keys=True, cls=CustomJSONEncoder)
                report_message = f"Schema migration results ({os.environ['DEPLOYMENT_STAGE']} env)"
                self._upload_to_slack("schema_migration_report.json", report_str, report_message)
                # Cleanup leftover schema migration files
                self.s3_provider.delete_prefix(artifact_bucket, self.get_key_prefix(f"schema_migration/{execution_id}"))

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
            response, _ = gather_collections()
        elif step_name == "collection_migrate":
            collection_id = os.environ["COLLECTION_ID"]
            collection_version_id = os.environ["COLLECTION_VERSION_ID"]
            collection_migrate = self.error_wrapper(self.collection_migrate, collection_id)
            response, _, _ = collection_migrate(
                collection_id=collection_id,
                collection_version_id=collection_version_id,
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
        elif step_name == "collection_cleanup":
            collection_version_id = os.environ["COLLECTION_VERSION_ID"]
            log_errors_and_cleanup = self.error_wrapper(self.log_errors_and_cleanup, collection_version_id)
            response = log_errors_and_cleanup(collection_version_id=collection_version_id)
        elif step_name == "report":
            response = self.report()
        self.logger.info("output", extra={"response": response})
        sfn_client = StepFunctionProvider().client
        sfn_client.send_task_success(
            taskToken=os.environ["TASK_TOKEN"], output=json.dumps(response, cls=CustomJSONEncoder)
        )
        return True
