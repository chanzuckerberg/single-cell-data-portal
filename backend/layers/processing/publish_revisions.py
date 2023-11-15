"""
After a migration has run, this can be used to publish all open revision for migrated collections.

$ aws batch submit-job --job-name publish-revisions \
  \ --job-queue <your_job_queue_ARN>
  \ --job-definition <your_job_definition_ARN>

"""
import logging
import os
from typing import Dict, List

from backend.layers.business.business import BusinessLogic
from backend.layers.business.entities import CollectionQueryFilter
from backend.layers.common.entities import (
    CollectionVersionWithDatasets,
    DatasetProcessingStatus,
)
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.processing.logger import configure_logging
from backend.layers.processing.process_logic import ProcessingLogic
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.schema_validator_provider import SchemaValidatorProvider
from backend.layers.thirdparty.uri_provider import UriProvider

configure_logging(level=logging.INFO)


class PublishRevisions(ProcessingLogic):
    def __init__(self, business_logic: BusinessLogic) -> None:
        super().__init__()
        self.business_logic = business_logic
        self.s3_provider = business_logic.s3_provider
        self.schema_validator = SchemaValidatorProvider()
        self.artifact_bucket = os.environ.get("ARTIFACT_BUCKET", "test-bucket")
        self.schema_version = self.schema_validator.get_current_schema_version()

    def check_datasets(self, collection_version: CollectionVersionWithDatasets) -> List[Dict]:
        """Check that all datasets have been migrated and are in a success state"""
        errors = []
        for dataset in collection_version.datasets:
            dataset_id = dataset.dataset_id.id
            dataset_version_id = dataset.version_id.id
            if not self.check_dataset_is_latest_schema_version(dataset):
                errors.append(
                    {
                        "message": "Dataset is not the latest schema version.",
                        "dataset_version_id": dataset_version_id,
                        "dataset_id": dataset_id,
                    }
                )
            elif dataset.status.processing_status != DatasetProcessingStatus.SUCCESS:
                errors.append(
                    {
                        "message": dataset.status.validation_message,
                        "dataset_status": dataset.status.to_dict(),
                        "dataset_version_id": dataset_version_id,
                        "dataset_id": dataset_id,
                    }
                )
        return errors

    def run(self):
        for collection_version in self.business_logic.get_collections(CollectionQueryFilter(is_published=False)):
            if collection_version.is_unpublished_version():
                _collection_version = self.business_logic.get_collection_version(collection_version.version_id)
                errors = self.check_datasets(_collection_version)
                if errors:
                    self.logger.error(
                        "Unable to publish collection version.",
                        extra={
                            "collection_id": collection_version.collection_id.id,
                            "collection_version_id": collection_version.version_id.id,
                            "errors": errors,
                        },
                    )
                else:
                    self.logger.info(
                        "Publishing collection version.",
                        extra={"collection_version_id": collection_version.version_id.id},
                    )
                    try:
                        self.business_logic.publish_collection_version(collection_version.version_id)
                    except Exception:
                        logging.exception("Failed to publish collection version.")


if __name__ == "__main__":
    business_logic = BusinessLogic(
        DatabaseProvider(),
        None,
        None,
        None,
        S3Provider(),
        UriProvider(),
    )

    PublishRevisions(business_logic).run()
