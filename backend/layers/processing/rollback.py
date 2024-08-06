"""
This batch job is meant for rollback of datasets or collections in response to migration failures--
either pre- or post- publish. It can either:
1) rollback all datasets to their previous dataset version in all private collections
2) rollback all datasets in an input list of private collections
3) rollback an input list of private datasets
4) rollback all public collections to their previous collection version
5) rollback an input list of public collections
depending on the input "ROLLBACK_TYPE" and optional "ENTITY_LIST" environment variables.

Example usages:

1) Rollback all datasets in all private collections:
$ aws batch submit-job --job-name rollback \
  --job-queue <your_job_queue_ARN> \
  --job-definition <your_job_definition_ARN> \
  --container-overrides '{
    "environment": [{"name": "ROLLBACK_TYPE", "value": "private_collections"}]
  }'

2) Rollback all datasets in an input list of private collections:
$ aws batch submit-job --job-name rollback \
  --job-queue <your_job_queue_ARN> \
  --job-definition <your_job_definition_ARN> \
  --container-overrides '{
    "environment": [{"name": "ROLLBACK_TYPE", "value": "private_collection_list"},
    {"name": "ENTITY_LIST", "value": "collection_version_id1,collection_version_id2"}]
  }'

3) Rollback an input list of private datasets:
$ aws batch submit-job --job-name rollback \
  --job-queue <your_job_queue_ARN> \
  --job-definition <your_job_definition_ARN> \
  --container-overrides '{
    "environment": [{"name": "ROLLBACK_TYPE", "value": "private_dataset_list"},
    {"name": "ENTITY_LIST", "value": "dataset_version_id1,dataset_version_id2"}]
  }'

4) Rollback all public collections to their previous collection version:
$ aws batch submit-job --job-name rollback \
  --job-queue <your_job_queue_ARN> \
  --job-definition <your_job_definition_ARN> \
  --container-overrides '{
    "environment": [{"name": "ROLLBACK_TYPE", "value": "public_collections"}]
  }'

5) Rollback an input list of public collections:
$ aws batch submit-job --job-name rollback \
  --job-queue <your_job_queue_ARN> \
  --job-definition <your_job_definition_ARN> \
  --container-overrides '{
    "environment": [{"name": "ROLLBACK_TYPE", "value": "public_collection_list"},
    {"name": "ENTITY_LIST", "value": "canonical_collection_id1,canonical_collection_id2"}]
  }'
"""

import os
from enum import Enum
from typing import List

from backend.layers.business.business import BusinessLogic, CollectionQueryFilter
from backend.layers.common.entities import (
    CollectionId,
    CollectionVersion,
    CollectionVersionId,
    DatasetVersion,
    DatasetVersionId,
)
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.uri_provider import UriProvider


class RollbackType(Enum):
    PRIVATE_COLLECTIONS = "private_collections"
    PUBLIC_COLLECTIONS = "public_collections"
    PUBLIC_COLLECTION_LIST = "public_collection_list"
    PRIVATE_COLLECTION_LIST = "private_collection_list"
    PRIVATE_DATASET_LIST = "private_dataset_list"


class RollbackEntity:
    def __init__(
        self, business_logic: BusinessLogic, rollback_type: RollbackType, entity_id_list: List[str] = None
    ) -> None:
        self.business_logic = business_logic
        self.rollback_type = rollback_type
        self.entity_id_list = entity_id_list

    def rollback(self):
        if self.rollback_type == RollbackType.PRIVATE_COLLECTIONS:
            self.rollback_private_collections()
        elif self.rollback_type == RollbackType.PUBLIC_COLLECTIONS:
            self.rollback_public_collections()
        elif self.rollback_type == RollbackType.PUBLIC_COLLECTION_LIST:
            self.rollback_public_collection_list()
        elif self.rollback_type == RollbackType.PRIVATE_COLLECTION_LIST:
            self.rollback_private_collection_list()
        elif self.rollback_type == RollbackType.PRIVATE_DATASET_LIST:
            self.rollback_private_dataset_list()
        else:
            raise ValueError(f"Invalid rollback type: {self.rollback_type}")

    def rollback_private_collections(self) -> None:
        """
        Rollback all the datasets in all private collections. This will restore the
        state of private collections to their pre-migration state. Then, triggers deletion of the DB
        references and S3 assets for the rolled back DatasetVersions.
        """
        filter = CollectionQueryFilter(is_published=False)
        collections = self.business_logic.get_collections(filter)
        rolled_back_datasets = []
        for collection in collections:
            collection_datasets = self.rollback_private_collection(collection.version_id)
            rolled_back_datasets.extend(collection_datasets)
        self._clean_up_rolled_back_datasets(rolled_back_datasets)

    def rollback_private_collection_list(self) -> None:
        """
        Rollback all the datasets from the input list of private collections. This will restore the
        state of the list of private collections to their pre-migration state. Then, triggers deletion of the DB
        references and S3 assets for the rolled back DatasetVersions.
        """
        rolled_back_datasets = []
        for collection_version_id in self.entity_id_list:
            collection_dataset_versions = self.rollback_private_collection(CollectionVersionId(collection_version_id))
            rolled_back_datasets.extend(collection_dataset_versions)
        self._clean_up_rolled_back_datasets(rolled_back_datasets)

    def rollback_private_collection(self, collection_version_id: CollectionVersionId) -> List[DatasetVersion]:
        """
        Rolls back the dataset versions for all datasets in the given private collection. This will restore the state of
        the private collection to its pre-migration state.
        """
        collection_version = self.business_logic.get_collection_version(collection_version_id)
        for dataset in collection_version.datasets:
            self.rollback_private_dataset(dataset.version_id, collection_version_id)
        return collection_version.datasets

    def rollback_private_dataset_list(self) -> None:
        """
        Rolls back the Datasets in the CollectionVersions associated with each DatasetVersionId passed in, to their
        respective previous, most recently created DatasetVersion. Then, triggers deletion of the DB references and S3
        assets for the rolled back DatasetVersions.
        """
        rolled_back_datasets = []
        for dataset_version_id in self.entity_id_list:
            dataset_version_id = DatasetVersionId(dataset_version_id)
            rolled_back_dataset = self.rollback_private_dataset(dataset_version_id)
            rolled_back_datasets.append(rolled_back_dataset)
        self._clean_up_rolled_back_datasets(rolled_back_datasets)

    def rollback_private_dataset(
        self, dataset_version_id: DatasetVersionId, collection_version_id: CollectionVersionId = None
    ) -> DatasetVersion:
        """
        For a given DatasetVersionId and unpublished CollectionVersionId, rolls back the associated Dataset in the
        CollectionVersion to its previous, most recently created DatasetVersion and deletes the given DatasetVersion.

        :param dataset_version_id: DatasetVersionId of the DatasetVersion to rollback
        :param collection_version_id: CollectionVersionId of the CollectionVersion to rollback the DatasetVersion
        in. If not passed in, the CollectionVersionId will be determined from the DatasetVersionId.
        :return: DatasetVersion that was rolled back
        """
        cv_id = collection_version_id
        dataset_version = self.business_logic.get_dataset_version(dataset_version_id)
        if cv_id is None:
            # account for collection potentially having a migration revision and a non-migration revision
            collection_versions = self.business_logic.get_unpublished_collection_versions_from_canonical(
                dataset_version.collection_id
            )
            for cv in collection_versions:
                if cv_id is not None:
                    break
                for dataset in cv.datasets:
                    if dataset.version_id == dataset_version_id:
                        cv_id = cv.version_id
                        break
        if cv_id is None:
            raise ValueError(f"An Associated CollectionVersion not found for DatasetVersion {dataset_version_id}")
        self.business_logic.restore_previous_dataset_version(cv_id, dataset_version.dataset_id)
        return dataset_version

    def rollback_public_collections(self) -> None:
        """
        Rollback each public collection to its previous, most recently published CollectionVersion (if one exists).
        This will restore the state of public collections to their pre-migration state. Then, triggers deletion of the DB
        references and S3 assets for the rolled back CollectionVersions and their associated DatasetVersions.
        """
        filter = CollectionQueryFilter(is_published=True)
        collections = self.business_logic.get_collections(filter)
        rolled_back_collection_versions = []
        for collection in collections:
            rolled_back_collection_version = self.rollback_public_collection(collection.collection_id)
            rolled_back_collection_versions.append(rolled_back_collection_version)
        self._clean_up_published_collection_versions(rolled_back_collection_versions)

    def rollback_public_collection_list(self) -> None:
        """
        Rollback each public collection in the input list to its previous, most recently published CollectionVersion (if
        one exists). This will restore the state of the list of public collections to their pre-migration state. Then,
        triggers deletion of the DB references and S3 assets for the rolled back CollectionVersions and their associated
        DatasetVersions.
        """
        rolled_back_collection_versions = []
        for collection_id in self.entity_id_list:
            rolled_back_collection_version = self.rollback_public_collection(CollectionId(collection_id))
            rolled_back_collection_versions.append(rolled_back_collection_version)
        self._clean_up_published_collection_versions(rolled_back_collection_versions)

    def rollback_public_collection(self, collection_id: CollectionId) -> CollectionVersion:
        """
        For a given public CollectionId, rolls back the associated Collection to its previous, most recently published
        CollectionVersion (if one exists). This will restore the state of the public collection to its pre-migration
        state.

        :param collection_id: CollectionId of the Collection to rollback
        :return: CollectionVersion that was rolled back
        """
        return self.business_logic.restore_previous_collection_version(collection_id)

    def _clean_up_rolled_back_datasets(self, rolled_back_datasets: List[DatasetVersion]) -> None:
        """
        Triggers deletion of the DB references and S3 assets for the rolled back DatasetVersions.
        """
        # TODO: replace this sync call with an async lambda call once available
        self.business_logic.delete_dataset_versions(rolled_back_datasets)

    def _clean_up_published_collection_versions(self, rolled_back_collection_versions: List[CollectionVersion]) -> None:
        """
        Triggers deletion of the DB references and S3 assets for the rolled back CollectionVersions' associated
        DatasetVersions, if they are not associated with any still-existing CollectionVersions.
        """
        datasets_to_rollback = []
        for rolled_back_collection_version in rolled_back_collection_versions:
            dataset_rollback_candidates = rolled_back_collection_version.datasets
            collection_version_history = self.business_logic.get_collection_versions_from_canonical(
                rolled_back_collection_version.collection_id
            )
            dataset_version_history = {dv.version_id.id: dv for cv in collection_version_history for dv in cv.datasets}
            for dataset_rollback_candidate in dataset_rollback_candidates:
                if dataset_rollback_candidate.version_id.id not in dataset_version_history:
                    datasets_to_rollback.append(dataset_rollback_candidate)
        self._clean_up_rolled_back_datasets(datasets_to_rollback)


if __name__ == "__main__":
    business_logic = BusinessLogic(
        DatabaseProvider(),
        None,
        None,
        None,
        S3Provider(),
        UriProvider(),
    )
    rollback_type_str = os.environ.get("ROLLBACK_TYPE", None)
    if rollback_type_str is None:
        raise ValueError("ROLLBACK_TYPE is required")

    rollback_type = RollbackType(rollback_type_str)
    if rollback_type in (
        RollbackType.PUBLIC_COLLECTION_LIST,
        RollbackType.PRIVATE_COLLECTION_LIST,
        RollbackType.PRIVATE_DATASET_LIST,
    ):
        entities_to_rollback = os.environ.get("ENTITY_LIST", None)
        entities_to_rollback_list = entities_to_rollback.split(",") if entities_to_rollback else None
        RollbackEntity(business_logic, rollback_type, entities_to_rollback_list).rollback()
    else:
        RollbackEntity(business_logic, rollback_type).rollback()
