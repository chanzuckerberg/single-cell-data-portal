import datetime
import uuid
from backend.corpora.common.entities import dataset
from backend.layers.persistence.persistence import DatabaseProviderInterface
from typing import Dict, List
from backend.layers.common.entities import (
    Collection,
    CollectionMetadata,
    CollectionVersion,
    Dataset,
    DatasetArtifact,
    DatasetMetadata,
    DatasetStatus,
    DatasetVersion,
)
import copy


class DatabaseProviderMock(DatabaseProviderInterface):

    """
    A mocked implementation for DatabaseProvider that uses in-memory dicts.
    This mock is to be used in tests only.
    Note: this implementation doesn't use immutability. Tests that assert immutability
    should NOT use this mock.
    """

    collections: Dict[str, str]
    collections_versions: Dict[str, CollectionVersion]
    datasets: Dict[str, str]
    datasets_versions: Dict[str, DatasetVersion]

    def __init__(self) -> None:
        super().__init__()
        self.collections = {}
        self.collections_versions = {}
        self.datasets = {}
        self.datasets_versions = {}

    @staticmethod
    def _id():
        return uuid.uuid4()

    def create_collection(self, collection_metadata: CollectionMetadata):
        collection_id = self._id()
        version_id = self._id()
        collection = Collection(collection_id, collection_metadata, None, [])
        version = CollectionVersion(version_id, self._id(), collection)
        self.collections_versions[version_id] = [version]
        self.collections[collection_id] = version_id

    def get_collection(self, collection_id: str) -> Collection:
        version = self.collections_versions[collection_id]
        return self.collections[version]

    def get_all_collections(self) -> List[Collection]:  # TODO: add filters if needed
        for version_id, collection in self.collections_versions:
            if version_id in self.collections.values():
                yield collection

    def delete_collection(self, collection_id: str) -> None:
        del self.collections[collection_id]

    def update_collection_metadata(self, version_id: str, collection_metadata: CollectionMetadata) -> None:
        self.collections_versions[version_id].collection.metadata = collection_metadata

    def update_collection_publisher_metadata(self, version_id: str, publisher_metadata: dict) -> None:
        self.collections_versions[version_id].collection.publisher_metadata = publisher_metadata

    def create_collection_version(self, collection_id: str) -> str:
        current_version = self.collections[collection_id]
        new_version_id = self._id()
        # Note: since datasets are immutable, there is no need to clone datasets here.
        new_version = copy.deepcopy(current_version)
        self.collections_versions[new_version_id] = CollectionVersion(new_version_id, new_version)
        return new_version_id

    def delete_collection_version(self, collection_id: str, version_id: str) -> None:
        del self.collections_versions[version_id]

    def get_collection_version(self, collection_id: str, version_id: str) -> CollectionVersion:
        return self.collections_versions[version_id]

    def publish_collection(self, collection_id: str, version_id: str) -> None:
        self.collections[collection_id] = version_id
        self.collections_versions[version_id].published_at = datetime.utcnow()

    def get_dataset(self, dataset_id: str) -> Dataset:
        version = self.datasets[dataset_id]
        return self.datasets_versions[version]

    def get_dataset_version(self, version_id: str) -> DatasetVersion:
        return self.datasets_versions[version_id]

    def get_all_datasets(self) -> List[Dataset]:  # TODO: add filters if needed
        for version_id, collection in self.datasets_versions:
            if version_id in self.collections.values():
                yield collection

    def get_dataset_artifacts(self, dataset_id: str) -> List[DatasetArtifact]:
        dataset = self.get_dataset(dataset_id)
        return dataset.artifacts

    def create_dataset(self, collection_version_id: str, dataset_metadata: DatasetMetadata) -> None:
        # Creates a dataset and initializes it with one version
        dataset_id = self._id()
        version_id = self._id()
        dataset = Dataset(dataset_id, None, dataset_metadata, [])
        version = DatasetVersion(version_id, dataset)
        self.dataset_versions[version_id] = [version]
        self.datasets[dataset_id] = version_id
        # Register the dataset to the original collection
        self.collections_versions[collection_version_id].collection.datasets.append(version_id)

    def create_dataset_version(self, dataset_id: str, dataset_metadata: DatasetMetadata):
        # Unused for now
        pass

    def create_dataset_artifact(self, version_id: str, artifact: DatasetArtifact) -> None:
        version = self.datasets_versions[version_id]
        version.dataset.artifacts.append(artifact)

    def update_dataset_status(self, version_id: str, status: DatasetStatus) -> None:
        dataset_version = self.dataset_versions[version_id]
        dataset_version.dataset.processing_status = status

    def add_dataset_to_collection_version(self, version_id: str, dataset_id: str) -> None:
        # Not needed for now - create_dataset does this
        # As an alternative, this could either be called by create_dataset
        pass

    def delete_dataset_from_collection_version(self, collection_version_id: str, dataset_version_id: str) -> None:
        version = self.collections_versions[collection_version_id]
        version.collection.datasets.remove(dataset_version_id)

    def replace_dataset_in_collection_version(
        self, collection_version_id: str, old_dataset_version_id: str, dataset_metadata: DatasetMetadata
    ) -> None:
        new_version_id = self._id()
        old_version = self.get_dataset_version(old_dataset_version_id)
        new_dataset = Dataset(old_version.dataset.id, None, dataset_metadata, [])
        new_version = DatasetVersion(new_version_id, new_dataset)
        self.datasets_versions[new_version_id] = new_version
        collection_version = self.collections_versions[collection_version_id]
        idx = collection_version.collection.datasets.index(old_dataset_version_id)
        collection_version.collection.datasets[idx] = new_version_id
