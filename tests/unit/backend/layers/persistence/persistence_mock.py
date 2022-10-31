from datetime import datetime
from importlib.metadata import metadata
from nntplib import ArticleInfo
import uuid
from backend.corpora.common.entities import dataset
from backend.layers.persistence.persistence import DatabaseProviderInterface
from typing import Dict, Iterable, List, Optional
from backend.layers.common.entities import (
    CollectionMetadata,
    CollectionVersion,
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
        self.collections = {} #rename to: active_collections
        self.collections_versions = {}
        self.datasets = {} # rename to: active_datasets
        self.datasets_versions = {}

    @staticmethod
    def _id():
        return str(uuid.uuid4())

    # TODO: add publisher_metadata here?
    def create_canonical_collection(self, owner: str, collection_metadata: CollectionMetadata) -> CollectionVersion:
        collection_id = self._id()
        version_id = self._id()
        version = CollectionVersion(collection_id, version_id, owner, collection_metadata, None, [], None)
        self.collections_versions[version_id] = version
        self.collections[collection_id] = version_id
        return version

    def get_collection(self, collection_id: str) -> CollectionVersion:
        version_id = self.collections[collection_id]
        return self.collections_versions[version_id]

    def get_all_collections_versions(self) -> Iterable[CollectionVersion]:  # TODO: add filters if needed
        for version_id, collection_version in self.collections_versions.items():
            if version_id in self.collections.values():
                yield collection_version

    def delete_collection(self, collection_id: str) -> None:
        del self.collections[collection_id]

    def save_collection_metadata(self, version_id: str, collection_metadata: CollectionMetadata) -> None:
        self.collections_versions[version_id].metadata = collection_metadata

    def save_collection_publisher_metadata(self, version_id: str, publisher_metadata: dict) -> None:
        self.collections_versions[version_id].publisher_metadata = publisher_metadata

    def add_collection_version(self, collection_id: str) -> str:
        current_version_id = self.collections[collection_id]
        current_version = self.collections_versions[current_version_id]
        new_version_id = self._id()
        # Note: since datasets are immutable, there is no need to clone datasets here.
        collection_version = CollectionVersion(
            collection_id=current_version.collection_id,
            version_id = new_version_id,
            owner = current_version.owner,
            metadata=current_version.metadata,
            publisher_metadata=current_version.publisher_metadata,
            datasets=current_version.datasets,
            published_at=None
        )
        self.collections_versions[new_version_id] = collection_version
        return new_version_id

    def delete_collection_version(self, version_id: str) -> None:
        # Can only delete an unpublished collection
        del self.collections_versions[version_id]

    def get_collection_version(self, version_id: str) -> CollectionVersion:
        return self.collections_versions[version_id]

    # MAYBE
    def finalize_collection_version(self, collection_id: str, version_id: str, published_at: Optional[datetime]) -> None:
        self.collections[collection_id] = version_id
        self.collections_versions[version_id].published_at = datetime.utcnow()

    # OR
    def update_collection_version_mapping(self, collection_id: str, version_id: str) -> None:
        self.collections[collection_id] = version_id

    def set_collection_version_published_at(self, version_id) -> None:
        self.collections_versions[version_id].published_at = datetime.utcnow()

    # END OR

    def get_dataset(self, dataset_id: str) -> DatasetVersion:
        version_id = self.datasets[dataset_id]
        return self.datasets_versions[version_id]

    def get_dataset_version(self, version_id: str) -> DatasetVersion:
        return self.datasets_versions[version_id]

    def get_all_datasets(self) -> Iterable[DatasetVersion]:  # TODO: add filters if needed
        for version_id, dataset_version in self.datasets_versions.items():
            if version_id in self.datasets.values():
                yield dataset_version

    def get_dataset_artifacts(self, dataset_id: str) -> List[DatasetArtifact]:
        dataset = self.get_dataset(dataset_id)
        return dataset.artifacts

    def create_canonical_dataset(self, collection_version_id: str, dataset_metadata: DatasetMetadata) -> DatasetVersion:
        # Creates a dataset and initializes it with one version
        dataset_id = self._id()
        version_id = self._id()
        version = DatasetVersion(
            dataset_id=dataset_id,
            version_id=version_id,
            processing_status=DatasetStatus.WAITING,
            metadata=dataset_metadata,
            artifacts=[]
        )
        self.datasets_versions[version_id] = version
        self.datasets[dataset_id] = version_id
        # Register the dataset to the original collection
        self.collections_versions[collection_version_id].datasets.append(version)
        return version

    def add_dataset_version(self, dataset_id: str, dataset_metadata: DatasetMetadata) -> str:
        # Unused for now
        pass

    def create_dataset_artifact(self, version_id: str, artifact: DatasetArtifact) -> None:
        version = self.datasets_versions[version_id]
        version.artifacts.append(artifact)

    def update_dataset_processing_status(self, version_id: str, status: DatasetStatus) -> None:
        dataset_version = self.datasets_versions[version_id]
        dataset_version.processing_status = status

    def add_dataset_to_collection_version(self, version_id: str, dataset_id: str) -> None:
        # Not needed for now - create_dataset does this
        # As an alternative, this could either be called by create_dataset
        pass

    def delete_dataset_from_collection_version(self, collection_version_id: str, dataset_version_id: str) -> None:
        version = self.collections_versions[collection_version_id]
        version.datasets = [d for d in version.datasets if d.version_id != dataset_version_id]

    def replace_dataset_in_collection_version(
        self, collection_version_id: str, old_dataset_version_id: str, dataset_metadata: DatasetMetadata
    ) -> None:
        new_version_id = self._id()
        old_version = self.get_dataset_version(old_dataset_version_id)
        new_version = DatasetVersion(
            dataset_id=old_version.dataset_id,
            version_id=new_version_id,
            processing_status=DatasetStatus.WAITING,
            metadata=dataset_metadata,
            artifacts=[]
        )
        self.datasets_versions[new_version_id] = new_version
        collection_version = self.collections_versions[collection_version_id]
        idx = next(i for i, e in enumerate(collection_version.datasets) if e.version_id == old_dataset_version_id)
        collection_version.datasets[idx] = new_version
