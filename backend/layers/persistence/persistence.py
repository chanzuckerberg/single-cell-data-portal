from datetime import date, datetime
from typing import List, Optional, Iterable
from backend.corpora.common.entities.dataset import Dataset
from backend.layers.common.entities import CollectionMetadata, CollectionVersion, DatasetArtifact, DatasetMetadata, DatasetStatus, DatasetVersion


class DatabaseProviderInterface:

    def initialize_canonical_collection(self, owner: str, collection_metadata: CollectionMetadata) -> CollectionVersion:
        """
        Initializes a canonical collection, generating a canonical collection_id and a new version_id.
        Returns the newly created CollectionVersion
        """
        pass


    def get_collection_version(self, version_id: str) -> CollectionVersion:
        """
        Retrieves a specific collection version by id
        """
        pass

    def get_all_version_for_collection(self, collection_id: str) -> List[CollectionVersion]:
        """
        Retrieves all versions for a specific collections, without filtering
        """
        pass

    def get_all_collections_versions(self) -> Iterable[CollectionVersion]:
        """
        Retrieves all versions of all collections.
        # TODO: add a filtering object.
        """
        pass

    def delete_canonical_collection(self, collection_id: str) -> None:
        """
        Deletes (tombstones) a canonical collection.
        """
        pass

    def save_collection_metadata(self, version_id: str, collection_metadata: CollectionMetadata) -> None:
        """
        Saves collection metadata for a collection version
        """
        pass

    def save_collection_publisher_metadata(self, version_id: str, publisher_metadata: dict) -> None:
        """
        Saves publisher metadata for a collection version
        """
        pass

    def add_collection_version(self, collection_id: str) -> CollectionVersion:
        """
        Adds a collection version to an existing canonical collection. The new version copies the following data from
         the previous version: owner, metadata, publisher_metadata, datasets (IDs).
        Returns the new version.
        """
        pass

    def delete_collection_version(self, version_id: str) -> None:
        """
        Deletes a collection version.
        """
        pass

    def finalize_collection_version(self, collection_id: str, version_id: str, published_at: Optional[datetime]) -> None:
        """
        Finalizes a collection version. This is equivalent to calling:
        1. update_collection_version_mapping
        2. set_collection_version_published_at
        """
        pass

    def update_collection_version_mapping(self, collection_id: str, version_id: str) -> None:
        """
        Updates the mapping between the canonical collection `collection_id` and its `version_id`
        """
        pass

    def set_collection_version_published_at(self, version_id: str, published_at: datetime) -> None:
        """
        Sets the `published_at` flag for a collection version
        """
        pass

    def get_dataset_version(self, dataset_version_id: str) -> DatasetVersion:
        """
        Returns a dataset version by id.
        """
        pass

    def get_all_versions_for_dataset(self, dataset_id: str) -> List[DatasetVersion]:
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

    def get_dataset_artifacts(self, dataset_version_id: str) -> List[DatasetArtifact]:
        """
        Returns all the artifacts for a specific dataset version
        """
        pass

    def initialize_canonical_dataset(self, collection_version_id: str, dataset_metadata: DatasetMetadata) -> DatasetVersion:
        """
        Initializes a canonical dataset, generating a dataset_id and a dataset_version_id.
        Returns the newly created DatasetVersion.
        """
        pass

    def add_dataset_version(self, dataset_id: str, dataset_metadata: DatasetMetadata) -> DatasetVersion:
        """
        Adds a new dataset version to a canonical dataset.
        Returns the newly created DatasetVersion.
        """
        pass

    def add_dataset_artifact(self, version_id: str, artifact: DatasetArtifact) -> None:
        """
        Adds a dataset artifact to an existing dataset version.
        """
        pass

    def update_dataset_processing_status(self, version_id: str, status: DatasetStatus) -> None:
        """
        Updates the processing status for a dataset version.
        """

    def add_dataset_to_collection_version_mapping(self, collection_version_id: str, dataset_version_id: str) -> None:
        """
        Adds a mapping between an existing collection version and a dataset version
        """
        pass

    def delete_dataset_from_collection_version(self, collection_version_id: str, dataset_version_id: str) -> None:
        """
        Removes a mapping between a collection version and a dataset version
        """
        pass

    def replace_dataset_in_collection_version(self, collection_version_id: str, old_dataset_version_id: str, dataset_metadata: DatasetMetadata) -> None:
        """
        Replaces an existing mapping between a collection version and a dataset version
        """
        pass
