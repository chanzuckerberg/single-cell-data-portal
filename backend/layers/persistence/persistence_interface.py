from datetime import datetime
from typing import Iterable, List, Optional

from backend.layers.common.entities import (
    CollectionId,
    CollectionMetadata,
    CollectionVersion,
    CollectionVersionId,
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
        pass

    def get_collection_version(self, version_id: CollectionVersionId) -> CollectionVersion:
        """
        Retrieves a specific collection version by id
        """
        pass

    def get_collection_mapped_version(self, collection_id: CollectionId) -> Optional[CollectionVersion]:
        """
        Retrieves the latest mapped version for a collection
        """
        pass

    def get_all_versions_for_collection(self, collection_id: CollectionId) -> List[CollectionVersion]:
        """
        Retrieves all versions for a specific collections, without filtering
        """
        pass

    def get_all_collections_versions(self) -> Iterable[CollectionVersion]:
        """
        Retrieves all versions of all collections.
        TODO: for performance reasons, it might be necessary to add a filtering parameter here.
        """
        pass

    def get_all_mapped_collection_versions(self) -> Iterable[CollectionVersion]:
        """
        Retrieves all the collection versions that are mapped to a canonical collection.
        """
        pass

    def delete_canonical_collection(self, collection_id: CollectionId) -> None:
        """
        Deletes (tombstones) a canonical collection.
        """
        pass

    def save_collection_metadata(
        self, version_id: CollectionVersionId, collection_metadata: CollectionMetadata
    ) -> None:
        """
        Saves collection metadata for a collection version
        """
        pass

    def save_collection_publisher_metadata(
        self, version_id: CollectionVersionId, publisher_metadata: Optional[dict]
    ) -> None:
        """
        Saves publisher metadata for a collection version. Specify None to remove it
        """
        pass

    def add_collection_version(self, collection_id: CollectionId) -> CollectionVersionId:
        """
        Adds a collection version to an existing canonical collection. The new version copies the following data from
         the previous version: owner, metadata, publisher_metadata, datasets (IDs).
        Returns the new version.
        """
        pass

    def delete_collection_version(self, version_id: CollectionVersionId) -> None:
        """
        Deletes a collection version.
        """
        pass

    def finalize_collection_version(
        self,
        collection_id: CollectionId,
        version_id: CollectionVersionId,
        published_at: Optional[datetime],
    ) -> None:
        """
        Finalizes a collection version. This is equivalent to calling:
        1. update_collection_version_mapping
        2. set_published_at
        """
        pass

    def update_collection_version_mapping(self, collection_id: CollectionId, version_id: CollectionVersionId) -> None:
        """
        Updates the mapping between the canonical collection `collection_id` and its `version_id`
        """
        pass

    def set_collection_version_published_at(self, version_id: CollectionVersionId, published_at: datetime) -> None:
        """
        Sets the `published_at` datetime for a collection version and its datasets
        """
        pass

    def get_dataset_version(self, dataset_version_id: DatasetVersionId) -> DatasetVersion:
        """
        Returns a dataset version by id.
        """
        pass

    def get_all_versions_for_dataset(self, dataset_id: DatasetId) -> List[DatasetVersion]:
        """
        Returns all dataset versions for a canonical dataset_id
        """
        pass

    def get_all_datasets(self) -> Iterable[DatasetVersion]:  # TODO: add filters if needed
        """
        Returns all dataset versions.
        # TODO: Add filtering
        """
        pass

    def get_dataset_artifacts_by_version_id(self, dataset_version_id: DatasetVersionId) -> List[DatasetArtifact]:
        """
        Returns all the artifacts for a specific dataset version
        """
        pass

    def create_canonical_dataset(self, collection_version_id: CollectionVersionId) -> DatasetVersion:
        """
        Initializes a canonical dataset, generating a dataset_id and a dataset_version_id.
        Returns the newly created DatasetVersion.
        """
        pass

    def add_dataset_artifact(
        self, version_id: DatasetVersionId, artifact_type: str, artifact_uri: str
    ) -> DatasetArtifactId:
        """
        Adds a dataset artifact to an existing dataset version.
        """
        pass

    def update_dataset_processing_status(self, version_id: DatasetVersionId, status: DatasetProcessingStatus) -> None:
        """
        Updates the processing status for a dataset version.
        """
        pass

    def update_dataset_validation_status(self, version_id: DatasetVersionId, status: DatasetValidationStatus) -> None:
        """
        Updates the validation status for a dataset version.
        """
        pass

    def update_dataset_upload_status(self, version_id: DatasetVersionId, status: DatasetUploadStatus) -> None:
        """
        Updates the upload status for a dataset version.
        """
        pass

    def update_dataset_conversion_status(
        self, version_id: DatasetVersionId, status_type: str, status: DatasetConversionStatus
    ) -> None:
        """
        Updates the conversion status for a dataset version and for `status_type`
        """
        pass

    def update_dataset_validation_message(self, version_id: DatasetVersionId, validation_message: str) -> None:
        """
        Updates the validation message for a dataset version
        """
        pass

    def get_dataset_version_status(self, version_id: DatasetVersionId) -> DatasetStatus:
        """
        Returns the status for a dataset version
        """
        pass

    def set_dataset_metadata(self, version_id: DatasetVersionId, metadata: DatasetMetadata) -> None:
        """
        Sets the metadata for a dataset version
        """
        pass

    def add_dataset_to_collection_version_mapping(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId
    ) -> None:
        """
        Adds a mapping between an existing collection version and a dataset version
        """
        pass

    def delete_dataset_from_collection_version(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId
    ) -> None:
        """
        Removes a mapping between a collection version and a dataset version
        """
        pass

    def replace_dataset_in_collection_version(
        self, collection_version_id: CollectionVersionId, old_dataset_version_id: DatasetVersionId
    ) -> DatasetVersion:
        """
        Replaces an existing mapping between a collection version and a dataset version
        """
        pass

    def get_dataset_mapped_version(self, dataset_id: DatasetId) -> Optional[DatasetVersion]:
        """
        Returns the dataset version mapped to a canonical dataset_id, or None if not existing
        """
        pass
