from dataclasses import dataclass
from datetime import datetime
import uuid
from backend.layers.persistence.persistence import DatabaseProviderInterface
from typing import Dict, Iterable, List, Optional
from backend.layers.common.entities import (
    CanonicalCollection,
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
import copy


@dataclass
class CanonicalCollectionPrivate:
    canonical_collection: CanonicalCollection
    mapped_version: CollectionVersionId

class DatabaseProviderMock(DatabaseProviderInterface):

    """
    A mocked implementation for DatabaseProvider that uses in-memory dicts.
    This mock is to be used in tests only.
    NOTE: this implementation doesn't use immutability. Tests that assert immutability
    should NOT use this mock.
    NOTE: this implementation uses copy.deepcopy for each returned entity. This is necessary
    since in-memory entities are just pointers, so updates will also modify the objects that were
    previously returned. In a database, this won't happen since changes are persisted to disk.
    Deep copying solves this and makes the mock behave like a real database.
    """

    # A mapping between canonical collection ids and collection versions.
    collections: Dict[str, CanonicalCollectionPrivate]

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
        canonical = CanonicalCollection(collection_id, None, False)
        version = CollectionVersion(collection_id, version_id, owner, collection_metadata, None, [], None, datetime.utcnow(), canonical)
        self.collections_versions[version_id.id] = version
        # Don't set mappings here - those will be set when publishing the collection!
        return copy.deepcopy(version)

    def get_collection_mapped_version(self, collection_id: CollectionId) -> Optional[CollectionVersion]:
        cc = self.collections.get(collection_id.id)
        if cc is not None:
            return copy.deepcopy(self.collections_versions[cc.mapped_version.id])

    def get_all_collections_versions(self) -> Iterable[CollectionVersion]:  # TODO: add filters if needed
        for version in self.collections_versions.values():
            yield copy.deepcopy(version)

    def get_all_mapped_collection_versions(self) -> Iterable[CollectionVersion]:  # TODO: add filters if needed
        for version_id, collection_version in self.collections_versions.items():
            if version_id in self.collections.values():
                yield copy.deepcopy(collection_version)

    def delete_collection(self, collection_id: CollectionId) -> None:
        del self.collections[collection_id.id]

    def save_collection_metadata(self, version_id: CollectionVersionId, collection_metadata: CollectionMetadata) -> None:
        self.collections_versions[version_id.id].metadata = copy.deepcopy(collection_metadata)

    def save_collection_publisher_metadata(self, version_id: CollectionVersionId, publisher_metadata: Optional[dict]) -> None:
        self.collections_versions[version_id.id].publisher_metadata = copy.deepcopy(publisher_metadata)

    def add_collection_version(self, collection_id: CollectionId) -> CollectionVersion:
        cc = self.collections[collection_id.id]
        current_version_id = cc.mapped_version
        current_version = self.collections_versions[current_version_id.id]
        new_version_id = CollectionVersionId(self._id())
        # Note: since datasets are immutable, there is no need to clone datasets here, 
        # but the list that contains datasets needs to be copied, since it's a pointer.
        new_dataset_list = copy.deepcopy(current_version.datasets)

        collection_version = CollectionVersion(
            collection_id=current_version.collection_id,
            version_id = new_version_id,
            owner = current_version.owner,
            metadata=current_version.metadata,
            publisher_metadata=current_version.publisher_metadata,
            datasets=new_dataset_list,
            published_at=None,
            created_at=datetime.utcnow(),
            canonical_collection=cc.canonical_collection
        )
        self.collections_versions[new_version_id.id] = collection_version
        return copy.deepcopy(collection_version)

    def delete_collection_version(self, version_id: CollectionVersionId) -> None:
        # Can only delete an unpublished collection
        del self.collections_versions[version_id.id]

    def get_collection_version(self, version_id: CollectionVersionId) -> CollectionVersion:
        return copy.deepcopy(self.collections_versions.get(version_id.id))

    def get_all_versions_for_collection(self, collection_id: CollectionId) -> Iterable[CollectionVersion]:
        # On a database, will require a secondary index on `collection_id` for an optimized lookup
        for collection_version in self.collections_versions.values():
            if collection_version.collection_id == collection_id:
                yield copy.deepcopy(collection_version)

    # MAYBE
    def finalize_collection_version(self, collection_id: CollectionId, version_id: CollectionVersionId, published_at: Optional[datetime]) -> None:
        cc = self.collections[collection_id.id]g
        cc.mapped_version = version_id
        now = datetime.utcnow()
        # If the canonical collection has never been published, set the field
        if cc.canonical_collection.published_at is None:
            cc.canonical_collection.published_at = now
        self.collections_versions[version_id.id].published_at = now

    # OR
    # def update_collection_version_mapping(self, collection_id: CollectionId, version_id: CollectionVersionId) -> None:
    #     self.collections[collection_id.id] = version_id.id

    # def set_collection_version_published_at(self, version_id: CollectionVersionId) -> None:
    #     self.collections_versions[version_id.id].published_at = datetime.utcnow()

    # END OR

    def get_dataset(self, dataset_id: DatasetId) -> DatasetVersion:
        version_id = self.datasets[dataset_id.id]
        return copy.deepcopy(self.datasets_versions[version_id])

    def get_dataset_version(self, version_id: DatasetVersionId) -> DatasetVersion:
        return copy.deepcopy(self.datasets_versions.get(version_id.id))

    def get_all_datasets(self) -> Iterable[DatasetVersion]:  
        """
        For now, this only returns all the active datasets, i.e. the datasets that belong to a published collection
        """
        active_collections = self.get_all_mapped_collection_versions()
        active_datasets = [i.version_id.id for s in [c.datasets for c in active_collections] for i in s]
        for version_id, dataset_version in self.datasets_versions.items():
            if version_id in active_datasets:
                yield copy.deepcopy(dataset_version)

    def _get_all_datasets(self) -> Iterable[DatasetVersion]:  
        """
        Returns all the mapped datasets. Currently unused
        """
        for version_id, dataset_version in self.datasets_versions.items():
            if version_id in self.datasets.values():
                yield copy.deepcopy(dataset_version)

    def get_dataset_artifacts(self, dataset_version_id: DatasetId) -> List[DatasetArtifact]:
        dataset = self.datasets_versions[dataset_version_id.id]
        return copy.deepcopy(dataset.artifacts)

    def create_canonical_dataset(self, collection_version_id: CollectionVersionId) -> DatasetVersion:
        # Creates a dataset and initializes it with one version
        dataset_id = DatasetId(self._id())
        version_id = DatasetVersionId(self._id())
        collection_version = self.collections_versions[collection_version_id.id]
        version = DatasetVersion(
            dataset_id=dataset_id,
            version_id=version_id,
            collection_id=collection_version.collection_id,
            status=DatasetStatus.empty(),
            metadata=None,
            artifacts=[],
            created_at=datetime.utcnow()
        )
        self.datasets_versions[version_id.id] = version
        self.datasets[dataset_id.id] = version_id.id
        # Register the dataset to the original collection
        collection_version.datasets.append(version)
        return copy.deepcopy(version)

    def add_dataset_version(self, dataset_id: DatasetId) -> str:
        # Unused for now
        raise NotImplementedError

    def add_dataset_artifact(self, version_id: DatasetVersionId, artifact_type: str, artifact_uri: str) -> DatasetArtifactId:
        version = self.datasets_versions[version_id.id]
        artifact_id = DatasetArtifactId(self._id())
        version.artifacts.append(DatasetArtifact(artifact_id, artifact_type, artifact_uri))
        return artifact_id

    def set_dataset_metadata(self, version_id: DatasetVersionId, metadata: DatasetMetadata) -> None:
        version = self.datasets_versions[version_id.id]
        version.metadata = copy.deepcopy(metadata)

    def update_dataset_processing_status(self, version_id: DatasetVersionId, status: DatasetProcessingStatus) -> None:
        dataset_version = self.datasets_versions[version_id.id]
        dataset_version.status.processing_status = copy.deepcopy(status)

    def update_dataset_validation_status(self, version_id: DatasetVersionId, status: DatasetValidationStatus) -> None:
        dataset_version = self.datasets_versions[version_id.id]
        dataset_version.status.validation_status = copy.deepcopy(status)

    def update_dataset_upload_status(self, version_id: DatasetVersionId, status: DatasetUploadStatus) -> None:
        dataset_version = self.datasets_versions[version_id.id]
        dataset_version.status.upload_status = copy.deepcopy(status)

    def update_dataset_conversion_status(self, version_id: DatasetVersionId, status_type: str, status: DatasetConversionStatus) -> None:
        dataset_version = self.datasets_versions[version_id.id]
        existing_status = dataset_version.status
        setattr(existing_status, status_type, copy.deepcopy(status))

    def add_dataset_to_collection_version(self, version_id: CollectionVersionId, dataset_id: DatasetId) -> None:
        # Not needed for now - create_dataset does this
        # As an alternative, this could either be called by create_dataset
        pass

    def delete_dataset_from_collection_version(self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId) -> None:
        version = self.collections_versions[collection_version_id.id]
        version.datasets = [d for d in version.datasets if d.version_id != dataset_version_id]

    def replace_dataset_in_collection_version(
        self, collection_version_id: CollectionVersionId, old_dataset_version_id: DatasetVersionId
    ) -> DatasetVersion:
        new_version_id = DatasetVersionId(self._id())
        old_version = self.get_dataset_version(old_dataset_version_id)
        collection_version = self.collections_versions[collection_version_id.id]
        new_version = DatasetVersion(
            dataset_id=old_version.dataset_id,
            version_id=new_version_id,
            collection_id=collection_version.collection_id,
            status=DatasetStatus.empty(),
            metadata=None,
            artifacts=[],
            created_at=datetime.utcnow()
        )
        self.datasets_versions[new_version_id.id] = new_version
        
        idx = next(i for i, e in enumerate(collection_version.datasets) if e.version_id == old_dataset_version_id)
        collection_version.datasets[idx] = new_version
        return copy.deepcopy(new_version)

    def get_dataset_version_status(self, version_id: DatasetVersionId) -> DatasetStatus:
        return copy.deepcopy(self.datasets_versions[version_id.id].status)

    def get_dataset_mapped_version(self, dataset_id: DatasetId) -> Optional[DatasetVersion]:
        version_id = self.datasets.get(dataset_id.id)
        if version_id is not None:
            return copy.deepcopy(self.datasets_versions[version_id])