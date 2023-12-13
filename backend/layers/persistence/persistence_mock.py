import copy
from datetime import datetime
from fnmatch import fnmatchcase
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

from backend.layers.business.exceptions import CollectionIsPublishedException, DatasetIsPublishedException
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
            schema_version=None,
            canonical_collection=canonical,
            datasets=[],
        )
        self.collections_versions[version_id.id] = version
        # Don't set mappings here - those will be set when publishing the collection!
        return copy.deepcopy(version)

    def _update_version_with_canonical(
        self,
        version: Union[CollectionVersion, CollectionVersionWithDatasets],
        update_datasets: bool = False,
        get_tombstoned: bool = False,
    ):
        """
        Private method that returns a version updated with the canonical collection.
        This is equivalent to a database double lookup (or join).
        Note that for methods that require a table scan, this should be optimized by holding the
        canonical_collections table in memory
        """
        copied_version = copy.deepcopy(version)
        if update_datasets:
            datasets_to_include = []
            for dataset_version_id in copied_version.datasets:
                if dataset_version := self.get_dataset_version(dataset_version_id, get_tombstoned=get_tombstoned):
                    dataset_version = self._update_dataset_version_with_canonical(dataset_version)
                    datasets_to_include.append(dataset_version)
            # Replace 'datasets' array of Dataset version ids with 'datasets' array of actual Dataset versions
            copied_version.datasets = datasets_to_include
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
        return copied_version

    def get_collection_mapped_version(self, collection_id: CollectionId) -> Optional[CollectionVersionWithDatasets]:
        cc = self.collections.get(collection_id.id)
        if cc is not None:
            return self.get_collection_version_with_datasets(cc.version_id, get_tombstoned=True)

    def get_all_collections_versions(
        self, get_tombstoned: bool = False
    ) -> Iterable[CollectionVersion]:  # TODO: add filters if needed
        for version in self.collections_versions.values():
            updated_version = self._update_version_with_canonical(version)
            if not get_tombstoned and updated_version.canonical_collection.tombstoned:
                continue
            included_dataset_version_ids = []
            for d_v_id in updated_version.datasets:
                if not get_tombstoned and self.datasets_versions[d_v_id.id].canonical_dataset.tombstoned:
                    continue
                included_dataset_version_ids.append(d_v_id)
            yield updated_version

    def get_canonical_collection(self, collection_id: CollectionId) -> Optional[CanonicalCollection]:
        return self.collections.get(collection_id.id, None)

    def get_all_mapped_collection_versions(
        self, get_tombstoned: bool = False
    ) -> Iterable[CollectionVersion]:  # TODO: add filters if needed
        for version_id, collection_version in self.collections_versions.items():
            if version_id in [c.version_id.id for c in self.collections.values()]:
                collection_id = collection_version.collection_id.id
                if not get_tombstoned and self.collections[collection_id].tombstoned:
                    continue
                yield self._update_version_with_canonical(collection_version)

    def tombstone_collection(self, collection_id: CollectionId) -> None:
        collection = self.collections[collection_id.id]
        collection.tombstoned = True
        # Tombstone Datasets individually as well
        for collection_version in self.collections_versions.values():
            if collection_version.collection_id == collection.id:
                for dataset_version in collection_version.datasets:
                    self.datasets[self.datasets_versions[dataset_version.id].dataset_id.id].tombstoned = True

    def resurrect_collection(self, collection_id: CollectionId, datasets_to_resurrect: Iterable[str]) -> None:
        """
        Untombstones a canonical collection and the explicitly-passed list of constituent Dataset ids. Constituent
        Datasets whose ids are not included in this list will remain tombstoned.
        """
        collection = self.collections[collection_id.id]
        collection.tombstoned = False
        # Untombstone Datasets individually as well
        for d_id in datasets_to_resurrect:
            self.datasets[d_id].tombstoned = False

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
            schema_version=None,
            canonical_collection=cc,
        )
        self.collections_versions[new_version_id.id] = collection_version
        return new_version_id

    def delete_unpublished_collection(self, collection_id: CollectionId) -> None:
        collection = self.collections.get(collection_id.id)
        if collection:
            if collection.originally_published_at is not None:
                raise CollectionIsPublishedException("Can only delete unpublished collections")
            del self.collections[collection_id.id]

    def delete_collection_version(self, version_id: CollectionVersionId) -> None:
        collection_version = self.collections_versions.get(version_id.id)
        if collection_version and collection_version.published_at is not None:
            raise CollectionIsPublishedException("Can only delete unpublished collections")
        del self.collections_versions[version_id.id]

    def delete_datasets(self, datasets: List[Union[DatasetId, CanonicalDataset]]) -> None:
        for d in datasets:
            d_id = d.id if isinstance(d, DatasetId) else d.dataset_id.id
            dataset = self.datasets.get(d_id)
            if dataset.published_at:
                raise DatasetIsPublishedException(f"Published Dataset {d_id} cannot be deleted")
            dataset_versions = list(
                filter(lambda dv: dv.dataset_id == dataset.dataset_id, self.datasets_versions.values())
            )
            self._delete_dataset_version_and_artifact_rows(dataset_versions, None)  # None in place of Session
            del self.datasets[d_id]

    def delete_dataset_versions(self, dataset_versions: List[Union[DatasetVersionId, DatasetVersion]]) -> None:
        ids = [d_v.id if isinstance(d_v, DatasetVersionId) else d_v.version_id.id for d_v in dataset_versions]
        dataset_version_rows = [self.datasets_versions.get(_id) for _id in ids]
        self._delete_dataset_version_and_artifact_rows(dataset_version_rows, None)  # None in place of Session

    def _delete_dataset_version_and_artifact_rows(
        self, dataset_version_rows: List[DatasetVersion], session: Any
    ) -> None:
        for d_v_row in dataset_version_rows:
            del self.datasets_versions[d_v_row.version_id.id]  # Artifacts live on DatasetVersion; they get deleted

    def get_collection_version(self, version_id: CollectionVersionId) -> CollectionVersion:
        version = self.collections_versions.get(version_id.id)
        if version is not None:
            return self._update_version_with_canonical(version)

    def get_all_versions_for_collection(
        self, collection_id: CollectionId, get_tombstoned: bool = False
    ) -> Iterable[CollectionVersionWithDatasets]:
        # On a database, will require a secondary index on `collection_id` for an optimized lookup
        versions = []
        for collection_version in self.collections_versions.values():
            if collection_version.collection_id == collection_id:
                versions.append(
                    self._update_version_with_canonical(
                        collection_version, update_datasets=True, get_tombstoned=get_tombstoned
                    )
                )
        return versions

    def get_collection_version_with_datasets(
        self, version_id: CollectionVersionId, get_tombstoned: bool = False
    ) -> CollectionVersionWithDatasets:
        version = self.collections_versions.get(version_id.id)
        if not version:
            return None
        version = self._update_version_with_canonical(version, update_datasets=True, get_tombstoned=get_tombstoned)
        if not get_tombstoned and version.canonical_collection.tombstoned:
            return None
        return version

    # MAYBE
    def finalize_collection_version(
        self,
        collection_id: CollectionId,
        version_id: CollectionVersionId,
        schema_version: str,
        published_at: Optional[datetime] = None,
        update_revised_at: bool = False,
    ) -> List[str]:

        published_at = published_at if published_at else datetime.utcnow()

        dataset_ids_for_new_collection_version = []
        version = self.collections_versions[version_id.id]
        for dataset_version_id in version.datasets:
            dataset_version = self.get_dataset_version(dataset_version_id)
            if self.datasets[dataset_version.dataset_id.id].published_at is None:
                self.datasets[dataset_version.dataset_id.id].published_at = published_at
            elif self.datasets[dataset_version.dataset_id.id].revised_at is None:
                self.datasets[dataset_version.dataset_id.id].revised_at = published_at
            dataset_version.canonical_dataset.dataset_version_id = dataset_version.version_id
            dataset_ids_for_new_collection_version.append(dataset_version.dataset_id.id)
        previous_collection = self.collections.get(collection_id.id)

        dataset_version_ids_to_delete_from_s3 = []
        if previous_collection is None:
            self.collections[collection_id.id] = CanonicalCollection(
                id=collection_id,
                version_id=version_id,
                originally_published_at=published_at,
                revised_at=None,
                tombstoned=False,
            )

        else:
            # Check to see if any Datasets are missing from new version, tombstone if so
            previous_dataset_version_ids = self.collections_versions[previous_collection.version_id.id].datasets
            previous_dataset_ids = [
                self.datasets_versions[d_v_id.id].dataset_id.id for d_v_id in previous_dataset_version_ids
            ]
            for previous_dataset_id in previous_dataset_ids:
                if previous_dataset_id not in dataset_ids_for_new_collection_version:
                    # Dataset has been removed and needs to be tombstoned
                    self.datasets[previous_dataset_id].tombstoned = True
                    for dataset_version in self.datasets_versions.values():
                        if dataset_version.dataset_id == previous_dataset_id:
                            dataset_version_ids_to_delete_from_s3.append(dataset_version.version_id.id)

            new_collection = copy.deepcopy(previous_collection)
            new_collection.version_id = version_id
            if update_revised_at:
                new_collection.revised_at = published_at
            self.collections[collection_id.id] = new_collection
        self.collections_versions[version_id.id].published_at = published_at
        self.collections_versions[version_id.id].schema_version = schema_version

        return dataset_version_ids_to_delete_from_s3

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

    def get_dataset_version(self, version_id: DatasetVersionId, get_tombstoned: bool = False) -> DatasetVersion:
        dataset_version = self.datasets_versions.get(version_id.id)
        if dataset_version:
            if not get_tombstoned and dataset_version.canonical_dataset.tombstoned:
                return None
            return self._update_dataset_version_with_canonical(dataset_version)

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

    def get_dataset_versions_by_id(
        self, ids: List[DatasetVersionId], get_tombstoned: bool = False
    ) -> List[DatasetVersion]:
        dataset_versions = []
        for dv_id in ids:
            dataset_version = self._update_dataset_version_with_canonical(self.datasets_versions[dv_id.id])
            if not get_tombstoned and dataset_version.canonical_dataset.tombstoned:
                continue
            dataset_versions.append(dataset_version)
        return dataset_versions

    def get_all_dataset_versions_for_collection(
        self, collection_id: CollectionId, from_date: datetime = datetime.min
    ) -> List[DatasetVersion]:
        """
        Get all Dataset versions -- published and unpublished -- for a canonical Collection
        """
        from_date = datetime.min if from_date is None else from_date
        dataset_versions = list(
            filter(
                lambda dv: dv.collection_id == collection_id and dv.created_at >= from_date,
                self.datasets_versions.values(),
            )
        )
        return [self._update_dataset_version_with_canonical(dv) for dv in dataset_versions]

    def get_most_recent_active_dataset_version(self, dataset_id: DatasetId) -> Optional[DatasetVersion]:
        """
        Returns the most recent, active Dataset version for a canonical dataset_id
        """
        dataset_versions = list(filter(lambda dv: dv.dataset_id == dataset_id, self.datasets_versions.values()))
        if not dataset_versions:
            return None
        dataset_versions_map = {str(dv.version_id): dv for dv in dataset_versions}
        collection_id = dataset_versions[0].collection_id
        cv_dataset_versions_ids = [
            cv.datasets
            for cv in sorted(
                filter(lambda cv: cv.collection_id == collection_id, self.collections_versions.values()),
                key=lambda cv: cv.created_at,
                reverse=True,
            )
        ]
        for cv_dataset_versions in cv_dataset_versions_ids:
            for dv_id in dataset_versions_map:
                if DatasetVersionId(dv_id) in cv_dataset_versions:
                    return self._update_dataset_version_with_canonical(dataset_versions_map.get(dv_id))

    def get_all_versions_for_dataset(self, dataset_id: DatasetId) -> List[DatasetVersion]:
        """
        Returns all dataset versions for a canonical dataset_id. ***AT PRESENT THIS FUNCTION IS NOT USED***
        """
        versions = []
        for dataset_version in self.datasets_versions.values():
            if dataset_version.dataset_id == dataset_id:
                versions.append(self._update_dataset_version_with_canonical(dataset_version))
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
        canonical_dataset = CanonicalDataset(dataset_id, None, False, None)
        version = DatasetVersion(
            dataset_id=dataset_id,
            version_id=version_id,
            collection_id=collection_version.collection_id,
            status=DatasetStatus.empty(),
            metadata=None,
            artifacts=[],
            created_at=datetime.utcnow(),
            canonical_dataset=canonical_dataset,
        )
        self.datasets_versions[version_id.id] = version
        self.datasets[dataset_id.id] = canonical_dataset
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

    def update_dataset_artifact(self, artifact_id: DatasetArtifactId, artifact_uri: str) -> None:
        found_artifact = False
        for version in self.datasets_versions.values():
            if found_artifact:
                break
            for artifact in version.artifacts:
                if artifact.id == artifact_id:
                    artifact.uri = artifact_uri
                    found_artifact = True
                    break

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
        self,
        collection_version_id: CollectionVersionId,
        old_dataset_version_id: DatasetVersionId,
        new_dataset_version_id: DatasetVersionId = None,
    ) -> DatasetVersion:
        old_version = self.get_dataset_version(old_dataset_version_id)
        collection_version = self.collections_versions[collection_version_id.id]
        if new_dataset_version_id is None:
            new_dataset_version_id = DatasetVersionId()
            new_dataset_version = DatasetVersion(
                dataset_id=old_version.dataset_id,
                version_id=new_dataset_version_id,
                collection_id=collection_version.collection_id,
                status=DatasetStatus.empty(),
                metadata=None,
                artifacts=[],
                created_at=datetime.utcnow(),
                canonical_dataset=old_version.canonical_dataset,
            )
            self.datasets_versions[new_dataset_version_id.id] = new_dataset_version
        else:
            new_dataset_version = self.get_dataset_version(new_dataset_version_id)
            if collection_version.collection_id != new_dataset_version.collection_id:
                raise ValueError(
                    f"Dataset version {new_dataset_version_id} does not belong to collection {collection_version.collection_id}"
                )

        idx = next(i for i, e in enumerate(collection_version.datasets) if e == old_dataset_version_id)
        collection_version.datasets[idx] = new_dataset_version_id
        return copy.deepcopy(new_dataset_version)

    def get_dataset_version_status(self, version_id: DatasetVersionId) -> DatasetStatus:
        return copy.deepcopy(self.datasets_versions[version_id.id].status)

    def get_dataset_mapped_version(
        self, dataset_id: DatasetId, get_tombstoned: bool = False
    ) -> Optional[DatasetVersion]:
        cd = self.datasets.get(dataset_id.id)
        if cd is not None:
            if not get_tombstoned and cd.tombstoned:
                return None
            version = None if cd.dataset_version_id is None else self.datasets_versions[cd.dataset_version_id.id]
            return None if version is None else self._update_dataset_version_with_canonical(version)

    def get_collection_versions_by_schema(self, schema_version: str, has_wildcards: bool) -> List[CollectionVersion]:
        if has_wildcards:
            schema_version = schema_version.replace("_", "?")
            collection_versions = [
                cv
                for cv in self.collections_versions.values()
                if cv.schema_version is not None and fnmatchcase(cv.schema_version, schema_version)
            ]
        else:
            collection_versions = [
                cv for cv in self.collections_versions.values() if cv.schema_version == schema_version
            ]
        return copy.deepcopy(collection_versions)

    def get_previous_dataset_version_id(self, dataset_id: DatasetId) -> Optional[DatasetVersionId]:
        """
        Returns the previously created dataset version for a dataset.
        """

        self.datasets.get(dataset_id.id)
        datasets = [d for d in self.datasets_versions.values() if d.dataset_id == dataset_id]
        datasets.sort(key=lambda d: d.created_at, reverse=True)
        try:
            previous_dataset = datasets[1]
        except IndexError:
            return None
        else:
            return copy.deepcopy(previous_dataset.version_id)
