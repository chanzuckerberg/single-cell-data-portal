from datetime import datetime
from typing import Iterable, List, Optional, Tuple

from backend.layers.common.entities import (
    CanonicalCollection,
    CanonicalDataset,
    CollectionId,
    CollectionMetadata,
    CollectionVersion,
    CollectionVersionId,
    CollectionVersionWithDatasets,
    DatasetArtifact,
    DatasetArtifactId,
    DatasetConversionStatus,
    DatasetId,
    DatasetMetadata,
    DatasetProcessingStatus,
    DatasetStatus,
    DatasetUploadStatus,
    DatasetValidationStatus,
    DatasetVersion,
    DatasetVersionId,
)


class PersistenceException(Exception):
    pass


class DatabaseProviderInterface:
    def create_canonical_collection(
        self, owner: str, curator_name: str, collection_metadata: CollectionMetadata
    ) -> CollectionVersion:
        """
        Creates a new canonical collection, generating a canonical collection_id and a new version_id.
        Returns the newly created CollectionVersion
        """

    def get_canonical_collection(self, collection_id: CollectionId) -> CanonicalCollection:
        """
        Return the canonical collection with id `collection_id`
        """

    def get_canonical_dataset(self, dataset_id: DatasetId) -> CanonicalDataset:
        """
        Return the canonical Dataset with id `dataset_id`
        """

    def get_collection_version(self, version_id: CollectionVersionId) -> CollectionVersion:
        """
        Retrieves a specific collection version by id
        """

    def get_collection_version_with_datasets(
        self, version_id: CollectionVersionId, get_tombstoned: bool
    ) -> CollectionVersionWithDatasets:
        """
        Retrieves a specific collection version by id, with datasets
        """

    def get_collection_mapped_version(self, collection_id: CollectionId) -> Optional[CollectionVersionWithDatasets]:
        """
        Retrieves the latest mapped version for a collection
        """

    def get_all_versions_for_collection(self, collection_id: CollectionId) -> List[CollectionVersionWithDatasets]:
        """
        Retrieves all versions for a specific collections, without filtering
        """

    def get_all_collections_versions(self, get_tombstoned: bool) -> Iterable[CollectionVersion]:
        """
        Retrieves all versions of all collections.
        TODO: for performance reasons, it might be necessary to add a filtering parameter here.
        """

    def get_all_mapped_collection_versions(self) -> Iterable[CollectionVersion]:
        """
        Retrieves all the collection versions that are mapped to a canonical collection.
        """

    def delete_canonical_collection(self, collection_id: CollectionId) -> None:
        """
        Deletes (tombstones) a canonical collection.
        """

    def save_collection_metadata(
        self, version_id: CollectionVersionId, collection_metadata: CollectionMetadata
    ) -> None:
        """
        Saves collection metadata for a collection version
        """

    def save_collection_publisher_metadata(
        self, version_id: CollectionVersionId, publisher_metadata: Optional[dict]
    ) -> None:
        """
        Saves publisher metadata for a collection version. Specify None to remove it
        """

    def add_collection_version(self, collection_id: CollectionId) -> CollectionVersionId:
        """
        Adds a collection version to an existing canonical collection. The new version copies the following data from
         the previous version: owner, metadata, publisher_metadata, datasets (IDs).
        Returns the new version.
        """

    def delete_collection_version(self, version_id: CollectionVersionId) -> None:
        """
        Deletes a collection version.
        """

    def finalize_collection_version(
        self,
        collection_id: CollectionId,
        version_id: CollectionVersionId,
        published_at: Optional[datetime] = None,
        update_revised_at: bool = False,
    ) -> List[str]:
        """
        Finalizes a collection version. This is equivalent to calling:
        1. update_collection_version_mapping
        2. set_published_at
        """

    def update_collection_version_mapping(self, collection_id: CollectionId, version_id: CollectionVersionId) -> None:
        """
        Updates the mapping between the canonical collection `collection_id` and its `version_id`
        """

    def set_collection_version_published_at(self, version_id: CollectionVersionId, published_at: datetime) -> None:
        """
        Sets the `published_at` datetime for a collection version and its datasets
        """

    def get_dataset_version(self, dataset_version_id: DatasetVersionId, get_tombstoned: bool) -> DatasetVersion:
        """
        Returns a dataset version by id.
        """

    def get_most_recent_active_dataset_version(self, dataset_id: DatasetId) -> Optional[DatasetVersion]:
        """
        Returns the most recent, active Dataset version for a canonical dataset_id
        """

    def get_all_versions_for_dataset(self, dataset_id: DatasetId) -> List[DatasetVersion]:
        """
        Returns all dataset versions for a canonical dataset_id
        """

    def get_all_mapped_datasets_and_collections(
        self,
    ) -> Tuple[List[DatasetVersion], List[CollectionVersion]]:  # TODO: add filters if needed
        """
        Returns all dataset versions.
        # TODO: Add filtering
        """

    def get_dataset_artifacts_by_version_id(self, dataset_version_id: DatasetVersionId) -> List[DatasetArtifact]:
        """
        Returns all the artifacts for a specific dataset version
        """

    def create_canonical_dataset(self, collection_version_id: CollectionVersionId) -> DatasetVersion:
        """
        Initializes a canonical dataset, generating a dataset_id and a dataset_version_id.
        Returns the newly created DatasetVersion.
        """

    def add_dataset_artifact(
        self, version_id: DatasetVersionId, artifact_type: str, artifact_uri: str
    ) -> DatasetArtifactId:
        """
        Adds a dataset artifact to an existing dataset version.
        """

    def update_dataset_processing_status(self, version_id: DatasetVersionId, status: DatasetProcessingStatus) -> None:
        """
        Updates the processing status for a dataset version.
        """

    def update_dataset_validation_status(self, version_id: DatasetVersionId, status: DatasetValidationStatus) -> None:
        """
        Updates the validation status for a dataset version.
        """

    def update_dataset_upload_status(self, version_id: DatasetVersionId, status: DatasetUploadStatus) -> None:
        """
        Updates the upload status for a dataset version.
        """

    def update_dataset_conversion_status(
        self, version_id: DatasetVersionId, status_type: str, status: DatasetConversionStatus
    ) -> None:
        """
        Updates the conversion status for a dataset version and for `status_type`
        """

    def update_dataset_validation_message(self, version_id: DatasetVersionId, validation_message: str) -> None:
        """
        Updates the validation message for a dataset version
        """

    def get_dataset_version_status(self, version_id: DatasetVersionId) -> DatasetStatus:
        """
        Returns the status for a dataset version
        """

    def set_dataset_metadata(self, version_id: DatasetVersionId, metadata: DatasetMetadata) -> None:
        """
        Sets the metadata for a dataset version
        """

    def add_dataset_to_collection_version_mapping(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId
    ) -> None:
        """
        Adds a mapping between an existing collection version and a dataset version
        """

    def delete_dataset_from_collection_version(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId
    ) -> None:
        """
        Removes a mapping between a collection version and a dataset version
        """

    def replace_dataset_in_collection_version(
        self,
        collection_version_id: CollectionVersionId,
        old_dataset_version_id: DatasetVersionId,
        new_dataset_version_id: DatasetVersionId = None,
    ) -> DatasetVersion:
        """
        Replaces an existing mapping between a collection version and a dataset version
        """

    def get_dataset_mapped_version(self, dataset_id: DatasetId, get_tombstoned: bool) -> Optional[DatasetVersion]:
        """
        Returns the dataset version mapped to a canonical dataset_id, or None if not existing
        """

    def get_collection_versions_by_schema(self, schema_version: str, has_wildcards: bool) -> List[CollectionVersion]:
        """
        Returns a list with all collection versions that match the given schema_version. schema_version may contain
         wildcards.
        """

    def get_previous_dataset_version_id(self, dataset_id: DatasetId) -> Optional[DatasetVersionId]:
        """
        Returns the previously created dataset version for a dataset.
        """
