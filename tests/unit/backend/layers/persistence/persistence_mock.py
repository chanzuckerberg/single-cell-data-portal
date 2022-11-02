from datetime import datetime
from importlib.metadata import metadata
from nntplib import ArticleInfo
import uuid
from backend.corpora.common.entities import dataset
from backend.layers.persistence.persistence import DatabaseProviderInterface
from typing import Dict, Iterable, List, Optional
from backend.layers.common.entities import (
    CollectionId,
    CollectionMetadata,
    CollectionVersion,
    CollectionVersionId,
    DatasetArtifact,
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
import copy


class DatabaseProviderMock(DatabaseProviderInterface):

    """
    A mocked implementation for DatabaseProvider that uses in-memory dicts.
    This mock is to be used in tests only.
    Note: this implementation doesn't use immutability. Tests that assert immutability
    should NOT use this mock.
    """

    # A mapping between canonical collection ids and collection versions.
    collections: Dict[str, str]

    # All the collection versions
    collections_versions: Dict[str, CollectionVersion]

    # A mapping between canonical dataset ids and dataset versions.
    datasets: Dict[str, str]

    # All the dataset versions
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
        collection_id = CollectionId(self._id())
        version_id = CollectionVersionId(self._id())
        version = CollectionVersion(collection_id, version_id, owner, collection_metadata, None, [], None)
        self.collections_versions[version_id.id] = version
        # Don't set mappings here - those will be set when publishing the collection!
        return version

    def get_collection_mapped_version(self, collection_id: CollectionId) -> Optional[CollectionVersion]:
        version_id = self.collections.get(collection_id.id)
        if version_id is not None:
            return self.collections_versions[version_id]

    def get_all_collections_versions(self) -> Iterable[CollectionVersion]:  # TODO: add filters if needed
        return self.collections_versions.values()

    def get_all_mapped_collection_versions(self) -> Iterable[CollectionVersion]:  # TODO: add filters if needed
        for version_id, collection_version in self.collections_versions.items():
            if version_id in self.collections.values():
                yield collection_version

    def delete_collection(self, collection_id: CollectionId) -> None:
        del self.collections[collection_id.id]

    def save_collection_metadata(self, version_id: CollectionVersionId, collection_metadata: CollectionMetadata) -> None:
        self.collections_versions[version_id.id].metadata = collection_metadata

    def save_collection_publisher_metadata(self, version_id: CollectionVersionId, publisher_metadata: Optional[dict]) -> None:
        self.collections_versions[version_id.id].publisher_metadata = publisher_metadata

    def add_collection_version(self, collection_id: CollectionId) -> str:
        current_version_id = self.collections[collection_id.id]
        current_version = self.collections_versions[current_version_id]
        new_version_id = CollectionVersionId(self._id())
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
        self.collections_versions[new_version_id.id] = collection_version
        return new_version_id.id

    def delete_collection_version(self, version_id: CollectionVersionId) -> None:
        # Can only delete an unpublished collection
        del self.collections_versions[version_id.id]

    def get_collection_version(self, version_id: CollectionVersionId) -> CollectionVersion:
        return self.collections_versions[version_id.id]

    # MAYBE
    def finalize_collection_version(self, collection_id: CollectionId, version_id: CollectionVersionId, published_at: Optional[datetime]) -> None:
        self.collections[collection_id.id] = version_id.id
        self.collections_versions[version_id.id].published_at = datetime.utcnow()

    # OR
    def update_collection_version_mapping(self, collection_id: CollectionId, version_id: CollectionVersionId) -> None:
        self.collections[collection_id.id] = version_id.id

    def set_collection_version_published_at(self, version_id: CollectionVersionId) -> None:
        self.collections_versions[version_id.id].published_at = datetime.utcnow()

    # END OR

    def get_dataset(self, dataset_id: DatasetId) -> DatasetVersion:
        version_id = self.datasets[dataset_id.id]
        return self.datasets_versions[version_id]

    def get_dataset_version(self, version_id: DatasetVersionId) -> DatasetVersion:
        return self.datasets_versions[version_id.id]

    def get_all_datasets(self) -> Iterable[DatasetVersion]:  # TODO: add filters if needed
        for version_id, dataset_version in self.datasets_versions.items():
            if version_id in self.datasets.values():
                yield dataset_version

    def get_dataset_artifacts(self, dataset_id: DatasetId) -> List[DatasetArtifact]:
        dataset = self.get_dataset(dataset_id)
        return dataset.artifacts

    def create_canonical_dataset(self, collection_version_id: CollectionVersionId) -> DatasetVersion:
        # Creates a dataset and initializes it with one version
        dataset_id = DatasetId(self._id())
        version_id = DatasetVersionId(self._id())
        version = DatasetVersion(
            dataset_id=dataset_id,
            version_id=version_id,
            status=DatasetStatus.empty(),
            metadata=None,
            artifacts=[]
        )
        self.datasets_versions[version_id.id] = version
        self.datasets[dataset_id.id] = version_id.id
        # Register the dataset to the original collection
        self.collections_versions[collection_version_id.id].datasets.append(version)
        return version

    def add_dataset_version(self, dataset_id: DatasetId) -> str:
        # Unused for now
        raise NotImplementedError

    def create_dataset_artifact(self, version_id: DatasetVersionId, artifact: DatasetArtifact) -> None:
        version = self.datasets_versions[version_id.id]
        version.artifacts.append(artifact)

    def update_dataset_processing_status(self, version_id: DatasetVersionId, status: DatasetProcessingStatus) -> None:
        dataset_version = self.datasets_versions[version_id.id]
        dataset_version.status.processing_status = status

    def update_dataset_validation_status(self, version_id: DatasetVersionId, status: DatasetValidationStatus) -> None:
        dataset_version = self.datasets_versions[version_id.id]
        dataset_version.status.validation_status = status

    def update_dataset_upload_status(self, version_id: DatasetVersionId, status: DatasetUploadStatus) -> None:
        dataset_version = self.datasets_versions[version_id.id]
        dataset_version.status.upload_status = status

    def update_dataset_conversion_status(self, version_id: DatasetVersionId, status_type: str, status: DatasetConversionStatus) -> None:
        dataset_version = self.datasets_versions[version_id.id]
        existing_status = dataset_version.status
        setattr(existing_status, status_type, status)

    def add_dataset_to_collection_version(self, version_id: CollectionVersionId, dataset_id: DatasetId) -> None:
        # Not needed for now - create_dataset does this
        # As an alternative, this could either be called by create_dataset
        pass

    def delete_dataset_from_collection_version(self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId) -> None:
        version = self.collections_versions[collection_version_id.id]
        version.datasets = [d for d in version.datasets if d.version_id != dataset_version_id]

    def replace_dataset_in_collection_version(
        self, collection_version_id: CollectionVersionId, old_dataset_version_id: DatasetVersionId
    ) -> None:
        new_version_id = DatasetVersionId(self._id())
        old_version = self.get_dataset_version(old_dataset_version_id)
        new_version = DatasetVersion(
            dataset_id=old_version.dataset_id,
            version_id=new_version_id,
            status=DatasetStatus.empty(),
            metadata=None,
            artifacts=[]
        )
        self.datasets_versions[new_version_id.id] = new_version
        collection_version = self.collections_versions[collection_version_id.id]
        idx = next(i for i, e in enumerate(collection_version.datasets) if e.version_id == old_dataset_version_id)
        collection_version.datasets[idx] = new_version
