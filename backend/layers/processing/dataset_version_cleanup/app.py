import logging
from typing import List, TypedDict

from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import (
    DatasetVersionId,
)
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.thirdparty.s3_provider import S3Provider

logger = logging.getLogger("processing")

DatasetVersionCleanupEvent = TypedDict("DatasetVersionCleanupEvent", {"dataset_version_ids": List[str]})


def dataset_version_cleanup_handler(event: DatasetVersionCleanupEvent, _context) -> None:
    """
    Lambda function invoked by the migration step function that deletes
    DatasetArtifacts + DatasetVersions for an input List of DataSetVersionIds
    :param event: Lambda's event object
    :param context: Lambda's context object
    :return:
    """
    business_logic = BusinessLogic(DatabaseProvider(), None, None, None, S3Provider(), None)

    dataset_version_ids = [DatasetVersionId(entity_id=id) for id in event.get("dataset_version_ids", [])]
    dataset_versions = business_logic.database_provider.get_dataset_versions_by_id(
        dataset_version_ids, get_tombstoned=False
    )
    business_logic.delete_dataset_versions(dataset_versions)
