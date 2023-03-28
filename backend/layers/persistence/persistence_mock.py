import copy
from datetime import datetime
from typing import Dict, Iterable, List, Optional, Tuple, Union

from backend.layers.business.exceptions import CollectionIsPublishedException
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
from backend.layers.persistence.persistence import DatabaseProviderInterface


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
    collections: Dict[str, CanonicalCollection]

    # All the collection versions
    collections_versions: Dict[str, CollectionVersion]

    # A mapping between canonical dataset ids and dataset versions.
    datasets: Dict[str, CanonicalDataset]

    # All the dataset versions
    datasets_versions: Dict[str, DatasetVersion]

    def __init__(self) -> None:
        super().__init__()
        self.collections = {}  # rename to: active_collections
        self.collections_versions = {}
        self.datasets = {}  # rename to: active_datasets
        self.datasets_versions = {}

    # TODO: add publisher_metadata here?
    def create_canonical_collection(
        self, owner: str, curator_name: str, collection_metadata: CollectionMetadata
    ) -> CollectionVersion:
        collection_id = CollectionId()
        version_id = CollectionVersionId()
        canonical = CanonicalCollection(collection_id, None, None, None, False)
        version = CollectionVersion(
            collection_id=collection_id,
            version_id=version_id,
            owner=owner,
            curator_name=curator_name,
            metadata=collection_metadata,
            publisher_metadata=None,
            published_at=None,
            created_at=datetime.utcnow(),
            canonical_collection=canonical,
            datasets=[],
        )
        self.collections_versions[version_id.id] = version
        # Don't set mappings here - those will be set when publishing the collection!
        return copy.deepcopy(version)

    def _update_version_with_canonical(
        self, version: Union[CollectionVersion, CollectionVersionWithDatasets], update_datasets: bool = False
    ):
        """
        Private method that returns a version updated with the canonical collection.
        This is equivalent to a database double lookup (or join).
        Note that for methods that require a table scan, this should be optimized by holding the
        canonical_collections table in memory
        """
        copied_version = copy.deepcopy(version)
        if update_datasets:
            copied_version.datasets = [self.get_dataset_version(dataset_id) for dataset_id in copied_version.datasets]
            # Hack for business logic that uses isinstance
            copied_version.__class__ = CollectionVersionWithDatasets
        cc = self.collections.get(version.collection_id.id)
        if cc is None:
            return copied_version
        copied_version.canonical_collection = cc
        return copied_version

    def _update_dataset_version_with_canonical(self, version: DatasetVersion):
        cd = self.datasets.get(version.dataset_id.id)
        if cd is None:
            return copy.deepcopy(version)
        copied_version = copy.deepcopy(version)
        copied_version.canonical_dataset = cd
        print("\n\nDATASET\n\n")
        print(copied_version.canonical_dataset)
        return copied_version

    def get_collection_mapped_version(self, collection_id: CollectionId) -> Optional[CollectionVersionWithDatasets]:
        cc = self.collections.get(collection_id.id)
        if cc is not None:
            return self.get_collection_version_with_datasets(cc.version_id)

    def get_all_collections_versions(
        self, get_tombstoned: bool = False
    ) -> Iterable[CollectionVersion]:  # TODO: add filters if needed
        for version in self.collections_versions.values():
            updated_version = self._update_version_with_canonical(version)
            if not get_tombstoned and updated_version.canonical_collection.tombstoned:
                continue
            yield updated_version

    def get_all_mapped_collection_versions(self) -> Iterable[CollectionVersion]:  # TODO: add filters if needed
        for version_id, collection_version in self.collections_versions.items():
            if version_id in [c.version_id.id for c in self.collections.values()]:
                collection_id = collection_version.collection_id.id
                if self.collections[collection_id].tombstoned:
                    continue
                yield self._update_version_with_canonical(collection_version)

    def delete_canonical_collection(self, collection_id: CollectionId) -> None:
        collection = self.collections[collection_id.id]
        collection.tombstoned = True

    def save_collection_metadata(
        self, version_id: CollectionVersionId, collection_metadata: CollectionMetadata
    ) -> None:
        self.collections_versions[version_id.id].metadata = copy.deepcopy(collection_metadata)

    def save_collection_publisher_metadata(
        self, version_id: CollectionVersionId, publisher_metadata: Optional[dict]
    ) -> None:
        self.collections_versions[version_id.id].publisher_metadata = copy.deepcopy(publisher_metadata)

    def add_collection_version(self, collection_id: CollectionId) -> CollectionVersionId:
        cc = self.collections[collection_id.id]
        current_version_id = cc.version_id
        current_version = self.collections_versions[current_version_id.id]
        new_version_id = CollectionVersionId()
        # Note: since datasets are immutable, there is no need to clone datasets here,
        # but the list that contains datasets needs to be copied, since it's a pointer.
        new_dataset_list = copy.deepcopy(current_version.datasets)

        collection_version = CollectionVersion(
            collection_id=current_version.collection_id,
            version_id=new_version_id,
            owner=current_version.owner,
            curator_name=current_version.curator_name,
            metadata=current_version.metadata,
            publisher_metadata=current_version.publisher_metadata,
            datasets=new_dataset_list,
            published_at=None,
            created_at=datetime.utcnow(),
            canonical_collection=cc,
        )
        self.collections_versions[new_version_id.id] = collection_version
        return new_version_id

    def delete_collection_version(self, version_id: CollectionVersionId) -> None:
        collection_version = self.collections_versions.get(version_id.id)
        if collection_version.published_at is not None:
            raise CollectionIsPublishedException("Can only delete unpublished collections")
        del self.collections_versions[version_id.id]

    def get_collection_version(self, version_id: CollectionVersionId) -> CollectionVersion:
        version = self.collections_versions.get(version_id.id)
        if version is not None:
            return self._update_version_with_canonical(version)

    def get_all_versions_for_collection(self, collection_id: CollectionId) -> Iterable[CollectionVersionWithDatasets]:
        # On a database, will require a secondary index on `collection_id` for an optimized lookup
        versions = []
        for collection_version in self.collections_versions.values():
            if collection_version.collection_id == collection_id:
                versions.append(self._update_version_with_canonical(collection_version, update_datasets=True))
        return versions

    def get_collection_version_with_datasets(self, version_id: CollectionVersionId) -> CollectionVersionWithDatasets:
        version = self.collections_versions.get(version_id.id)
        if version is not None:
            return self._update_version_with_canonical(version, update_datasets=True)

    # MAYBE
    def finalize_collection_version(
        self,
        collection_id: CollectionId,
        version_id: CollectionVersionId,
        published_at: Optional[datetime] = None,
        update_revised_at: bool = False,
    ) -> None:

        published_at = published_at if published_at else datetime.utcnow()

        version = self.collections_versions[version_id.id]
        for dataset_version_id in version.datasets:
            dataset_version = self.get_dataset_version(dataset_version_id)
            if self.datasets[dataset_version.dataset_id.id].published_at is None:
                self.datasets[dataset_version.dataset_id.id].published_at = published_at
            if self.datasets[dataset_version.dataset_id.id].revised_at is None:
                self.datasets[dataset_version.dataset_id.id].revised_at = published_at
            dataset_version.canonical_dataset.dataset_version_id = dataset_version.version_id
        cc = self.collections.get(collection_id.id)
        if cc is None:
            self.collections[collection_id.id] = CanonicalCollection(
                id=collection_id,
                version_id=version_id,
                originally_published_at=published_at,
                revised_at=None,
                tombstoned=False,
            )

        else:
            new_cc = copy.deepcopy(cc)
            new_cc.version_id = version_id
            if update_revised_at:
                new_cc.revised_at = published_at
            self.collections[collection_id.id] = new_cc
        self.collections_versions[version_id.id].published_at = published_at

    # OR
    # def update_collection_version_mapping(self, collection_id: CollectionId, version_id: CollectionVersionId) -> None:
    #     self.collections[collection_id.id] = version_id.id

    # def set_collection_version_published_at(self, version_id: CollectionVersionId) -> None:
    #     self.collections_versions[version_id.id].published_at = datetime.utcnow()

    # END OR

    def get_canonical_dataset(self, dataset_id: DatasetId) -> CanonicalDataset:
        if dataset_id.id is None or dataset_id.id not in self.datasets:
            return None
        return self.datasets[dataset_id.id]

    def get_dataset_version(self, version_id: DatasetVersionId) -> DatasetVersion:
        version = self.datasets_versions.get(version_id.id)
        if version is not None:
            return self._update_dataset_version_with_canonical(version)

    def get_all_mapped_datasets_and_collections(self) -> Tuple[List[DatasetVersion], List[CollectionVersion]]:
        """
        For now, this only returns all the active datasets, i.e. the datasets that belong to a published collection
        """
        active_collections = list(self.get_all_mapped_collection_versions())
        active_datasets_ids = [i.id for s in [c.datasets for c in active_collections] for i in s]
        active_datasets = []
        for version_id, dataset_version in self.datasets_versions.items():
            if version_id in active_datasets_ids:
                active_datasets.append(self._update_dataset_version_with_canonical(dataset_version))
        return active_datasets, active_collections

    def _get_datasets(self, ids: List[DatasetVersionId]) -> List[DatasetVersion]:
        dataset_versions = []
        for dv_id in ids:
            dataset_versions.append(self._update_dataset_version_with_canonical(self.datasets_versions[dv_id.id]))
        return dataset_versions

    def get_all_mapped_collection_versions_with_datasets(self) -> Tuple[List[DatasetVersion], List[CollectionVersion]]:
        """
        Returns all mapped collection versions with their datasets
        """
        return self.get_all_mapped_datasets_and_collections()

    def get_all_versions_for_dataset(self, dataset_id: DatasetId) -> List[DatasetVersion]:
        versions = []
        for dataset_version in self.datasets_versions.values():
            if dataset_version.dataset_id == dataset_id:
                versions.append(dataset_version)
        return versions

    def _get_all_datasets(self) -> Iterable[DatasetVersion]:
        """
        Returns all the mapped datasets. Currently unused
        """
        for version_id, dataset_version in self.datasets_versions.items():
            if version_id in self.datasets.values():
                yield self._update_dataset_version_with_canonical(dataset_version)

    def get_dataset_artifacts_by_version_id(self, dataset_version_id: DatasetId) -> List[DatasetArtifact]:
        dataset = self.datasets_versions[dataset_version_id.id]
        return copy.deepcopy(dataset.artifacts)

    def create_canonical_dataset(self, collection_version_id: CollectionVersionId) -> DatasetVersion:
        # Creates a dataset and initializes it with one version
        dataset_id = DatasetId()
        version_id = DatasetVersionId()
        collection_version = self.collections_versions[collection_version_id.id]
        version = DatasetVersion(
            dataset_id=dataset_id,
            version_id=version_id,
            collection_id=collection_version.collection_id,
            status=DatasetStatus.empty(),
            metadata=None,
            artifacts=[],
            created_at=datetime.utcnow(),
            canonical_dataset=CanonicalDataset(dataset_id, None, None),
        )
        self.datasets_versions[version_id.id] = version
        self.datasets[dataset_id.id] = CanonicalDataset(
            dataset_id=dataset_id, dataset_version_id=None, published_at=None
        )
        return copy.deepcopy(version)

    def add_dataset_to_collection_version_mapping(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId
    ) -> None:
        self.collections_versions[collection_version_id.id].datasets.append(dataset_version_id)

    def add_dataset_artifact(
        self, version_id: DatasetVersionId, artifact_type: str, artifact_uri: str
    ) -> DatasetArtifactId:
        version = self.datasets_versions[version_id.id]
        artifact_id = DatasetArtifactId()
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

    def update_dataset_conversion_status(
        self, version_id: DatasetVersionId, status_type: str, status: DatasetConversionStatus
    ) -> None:
        dataset_version = self.datasets_versions[version_id.id]
        existing_status = dataset_version.status
        setattr(existing_status, status_type, copy.deepcopy(status))

    def update_dataset_validation_message(self, version_id: DatasetVersionId, validation_message: str) -> None:
        dataset_version = self.datasets_versions[version_id.id]
        dataset_version.status.validation_message = validation_message

    def add_dataset_to_collection_version(self, version_id: CollectionVersionId, dataset_id: DatasetId) -> None:
        # Not needed for now - create_dataset does this
        # As an alternative, this could either be called by create_dataset
        pass

    def delete_dataset_from_collection_version(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId
    ) -> None:
        version = self.collections_versions[collection_version_id.id]
        version.datasets = [d for d in version.datasets if d != dataset_version_id]

    def replace_dataset_in_collection_version(
        self, collection_version_id: CollectionVersionId, old_dataset_version_id: DatasetVersionId
    ) -> DatasetVersion:
        new_version_id = DatasetVersionId()
        old_version = self.get_dataset_version(old_dataset_version_id)
        collection_version = self.collections_versions[collection_version_id.id]
        new_version = DatasetVersion(
            dataset_id=old_version.dataset_id,
            version_id=new_version_id,
            collection_id=collection_version.collection_id,
            status=DatasetStatus.empty(),
            metadata=None,
            artifacts=[],
            created_at=datetime.utcnow(),
            canonical_dataset=old_version.canonical_dataset,
        )
        self.datasets_versions[new_version_id.id] = new_version

        idx = next(i for i, e in enumerate(collection_version.datasets) if e == old_dataset_version_id)
        collection_version.datasets[idx] = new_version_id
        return copy.deepcopy(new_version)

    def get_dataset_version_status(self, version_id: DatasetVersionId) -> DatasetStatus:
        return copy.deepcopy(self.datasets_versions[version_id.id].status)

    def get_dataset_mapped_version(self, dataset_id: DatasetId) -> Optional[DatasetVersion]:
        cd = self.datasets.get(dataset_id.id)
        if cd is not None:
            version = None if cd.dataset_version_id is None else self.datasets_versions[cd.dataset_version_id.id]
            return None if version is None else self._update_dataset_version_with_canonical(version)
