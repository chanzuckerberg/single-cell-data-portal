from typing import List
from backend.corpora.common.entities.dataset import Dataset
from backend.layers.common.entities import Collection, CollectionMetadata, CollectionVersion, DatasetAsset, DatasetMetadata, DatasetStatus


class DatabaseProviderInterface:
    def create_collection(self, collection_metadata: CollectionMetadata):
        pass

    def get_collection(self, collection_id: str) -> Collection:
        pass

    def get_all_collections(self) -> List[Collection]:  # TODO: add filters if needed
        pass

    def delete_collection(self, collection_id: str) -> None:
        pass

    def update_collection_metadata(self, version_id: str, collection_metadata: CollectionMetadata) -> None:
        pass

    def update_collection_publisher_metadata(self, version_id: str, publisher_metadata: dict) -> None:
        pass

    def create_collection_version(self, collection_id: str) -> str:
        pass

    def delete_collection_version(self, collection_id: str, version_id: str) -> None:
        pass

    def get_collection_version(self, collection_id: str, version_id: str) -> CollectionVersion:
        pass

    def publish_collection(self, collection_id: str, version_id: str) -> None:
        pass

    def get_dataset(self, dataset_id: str) -> Dataset:
        pass

    def get_dataset_version(self, version_id: str) -> DatasetVersion:
        pass

    def get_all_datasets(self) -> List[Dataset]:  # TODO: add filters if needed
        pass

    def get_dataset_artifacts(self, dataset_id: str) -> List[DatasetArtifact]:
        pass

    def create_dataset(self, collection_version_id: str, dataset_metadata: DatasetMetadata) -> None:
        pass

    def create_dataset_version(self, dataset_id: str, dataset_metadata: DatasetMetadata):
        pass

    def create_dataset_artifact(self, version_id: str, artifact: DatasetArtifact) -> None:
        pass

    def update_dataset_status(self, version_id: str, status: DatasetStatus) -> None:
        pass

    def add_dataset_to_collection_version(self, version_id: str, dataset_id: str) -> None:
        pass

    def delete_dataset_from_collection_version(self, collection_version_id: str, dataset_version_id: str) -> None:
        pass

    def replace_dataset_in_collection_version(