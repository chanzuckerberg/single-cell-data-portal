from datetime import datetime
from typing import List
from backend.corpora.common.entities.dataset import Dataset
from backend.layers.common.entities import CollectionMetadata, CollectionVersion, DatasetArtifact, DatasetMetadata, DatasetStatus


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

    def get_all_collections(self) -> Iterable[CollectionVersion]:
        """
        Retrieves all collections.
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
        Adds a collection version to an existing canonical collection.
        Returns the new version.
        """
        pass

    def delete_collection_version(self, version_id: str) -> None:
        """
        Deletes a collection version.
        """
        pass

    def finalize_collection_version(self, collection_id: str, version_id: str) -> None:
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

    def get_dataset(self, dataset_id: str) -> DatasetVersion:
        pass

    def get_dataset_version(self, version_id: str) -> DatasetVersion:
        pass

    def get_all_datasets(self) -> Iterable[DatasetVersion]:  # TODO: add filters if needed
        pass

    def get_dataset_artifacts(self, dataset_id: str) -> List[DatasetArtifact]:
        pass

    def create_dataset(self, collection_version_id: str, dataset_metadata: DatasetMetadata) -> None:
        pass

    def create_dataset_version(self, dataset_id: str, dataset_metadata: DatasetMetadata) -> str:
        pass

    def create_dataset_artifact(self, version_id: str, artifact: DatasetArtifact) -> None:
        pass

    def update_dataset_processing_status(self, version_id: str, status: DatasetStatus) -> None:
        pass

    def add_dataset_to_collection_version(self, version_id: str, dataset_id: str) -> None:
        pass

    def delete_dataset_from_collection_version(self, collection_version_id: str, dataset_version_id: str) -> None:
        pass

    def replace_dataset_in_collection_version(self, collection_version_id: str, old_dataset_version_id: str, dataset_metadata: DatasetMetadata) -> None:
        pass
