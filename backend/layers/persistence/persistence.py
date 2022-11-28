from datetime import datetime
from typing import List, Optional, Iterable
import uuid
from backend.layers.common.entities import (
    CanonicalCollection,
    CanonicalDataset,
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
from backend.layers.persistence.orm import (
    Collection as CollectionRow,
    CollectionVersion as CollectionVersionRow,
    Dataset as DatasetRow,
    DatasetVersion as DatasetVersionRow,
    DatasetArtifact as DatasetArtifactRow
)
from backend.layers.persistence.db_session import db_session_manager, _db_session_maker
from sqlalchemy import select

from backend.layers.persistence.persistence_interface import DatabaseProviderInterface


class DatabaseProvider(DatabaseProviderInterface):

    def __init__(self, **kwargs) -> None:
        self.db_session_manager = db_session_manager(_db_session_maker, **kwargs)

    @staticmethod
    def _generate_id():
        return str(uuid.uuid4())

    def _hydrate_dataset_version(self, dataset_version: DatasetVersionRow) -> DatasetVersion:
        """
        Populates canonical_dataset, artifacts, and status for DatasetVersionRow
        """
        if not dataset_version.canonical_dataset:
            dataset_version.canonical_dataset = self.get_canonical_dataset(dataset_version.dataset_id)
        if not dataset_version.artifacts:
            dataset_version.artifacts = self.get_dataset_artifacts(dataset_version.artifacts)
        return dataset_version

    def get_canonical_collection(self, collection_id: CollectionId) -> CanonicalCollection:
        with self.db_session_manager as session:
            collection = session.query(CollectionRow).filter_by(id=collection_id).one()
        return collection

    def get_canonical_dataset(self, dataset_id: DatasetId) -> CanonicalDataset:
        with self.db_session_manager as session:
            dataset = session.query(DatasetRow).filter_by(dataset_id=dataset_id).one()
        return dataset

    def create_canonical_collection(self, owner: str, collection_metadata: CollectionMetadata) -> CollectionVersion:
        """
        Creates a new canonical collection, generating a canonical collection_id and a new version_id.
        Returns the newly created CollectionVersion
        """
        collection_id = CollectionId((self._generate_id()))
        version_id = CollectionVersionId((self._generate_id()))
        canonical_collection = CollectionRow(id=collection_id,
                                             version_id=None,
                                             tombstoned=False,
                                             originally_published_at=None,
                                             )
        collection_version_row = CollectionVersionRow(collection_id=collection_id,
                                                      version_id=version_id,
                                                      owner=owner,
                                                      metadata=collection_metadata,
                                                      publisher_metadata=None,
                                                      published_at=None,
                                                      datasets=list(),
                                                      canonical_collection=canonical_collection,
                                                      created_at=datetime.utcnow(),
                                                      )
        with self.db_session_manager as session:
            session.add(canonical_collection)
            session.add(collection_version_row)

        return collection_version_row

    def get_collection_version(self, version_id: CollectionVersionId) -> CollectionVersion:
        """
        Retrieves a specific collection version by id
        """
        with self.db_session_manager as session:
            collection_version = session.query(CollectionVersionRow).filter_by(version_id=version_id).one()
        collection_version.canonical_collection = self.get_canonical_collection(collection_version.collection_id)
        return collection_version

    def get_collection_mapped_version(self, collection_id: CollectionId) -> Optional[CollectionVersion]:
        """
        Retrieves the latest mapped version for a collection
        """
        with self.db_session_manager as session:
            version_id = session.query(CollectionRow.version_id).filter_by(collection_id=collection_id).one_or_none() # noqa
            collection_version = session.query(CollectionVersionRow).filter_by(version_id=version_id).one_or_none()
        collection_version.canonical_collection = self.get_canonical_collection(collection_id)
        return collection_version

    def get_all_versions_for_collection(self, collection_id: CollectionId) -> List[CollectionVersion]:
        """
        Retrieves all versions for a specific collections, without filtering
        """
        with self.db_session_manager as session:
            versions = session.query(CollectionVersionRow).filter_by(collection_id=collection_id).all()
        canonical_collection = self.get_canonical_collection(collection_id)
        for i in range(len(versions)):
            versions[i].canonical_collection = canonical_collection
        return versions

    def get_all_collections_versions(self, get_tombstoned: bool = False) -> Iterable[CollectionVersion]:
        """
        Retrieves all versions of all collections.
        TODO: for performance reasons, it might be necessary to add a filtering parameter here.
        """
        with self.db_session_manager as session:
            if get_tombstoned:
                versions = session.query(CollectionVersionRow).all()
            else:
                version_ids = session.query(CollectionRow.version_id).filter_by(tombstoned=False).all()
                versions = session.query(CollectionVersionRow).filter(CollectionVersionRow.version_id in version_ids).all() # noqa
        # TODO: do we need to hydrate versions with canonical collections? would require a join or many lookup calls
        return versions

    def get_all_mapped_collection_versions(self, get_tombstoned: bool = False) -> Iterable[CollectionVersion]:
        """
        Retrieves all the collection versions that are mapped to a canonical collection.
        """
        with self.db_session_manager as session:
            if get_tombstoned:
                mapped_version_ids = session.query(CollectionRow.version_id).filter(CollectionRow.version_id.isnot(None)).all() # noqa
            else:
                mapped_version_ids = session.query(CollectionRow.version_id)\
                    .filter(CollectionRow.version_id.isnot(None))\
                    .filter_by(tombstoned=False)\
                    .all()

            versions = session.query(CollectionVersionRow).filter(CollectionVersionRow.version_id in mapped_version_ids).all() # noqa
        # TODO: do we need to hydrate versions with canonical collections? would require a join or many lookup calls
        return versions

    def delete_canonical_collection(self, collection_id: CollectionId) -> None:
        """
        Deletes (tombstones) a canonical collection.
        """
        with self.db_session_manager as session:
            canonical_collection = session.query(CollectionRow).filter_by(id=collection_id).one()
            canonical_collection.tombstoned = True

    def save_collection_metadata(
        self, version_id: CollectionVersionId, collection_metadata: CollectionMetadata
    ) -> None:
        """
        Saves collection metadata for a collection version
        """
        with self.db_session_manager as session:
            version = session.query(CollectionVersionRow).filter_by(version_id=version_id).one()
            version.metadata = collection_metadata

    def save_collection_publisher_metadata(
            self, version_id: CollectionVersionId, publisher_metadata: Optional[dict]
    ) -> None:
        """
        Saves publisher metadata for a collection version. Specify None to remove it
        """
        with self.db_session_manager as session:
            version = session.query(CollectionVersionRow).filter_by(verison_id=version_id)
            version.publisher_metadata = publisher_metadata

    def add_collection_version(self, collection_id: CollectionId) -> CollectionVersion:
        """
        Adds a collection version to an existing canonical collection. The new version copies the following data from
         the previous version: owner, metadata, publisher_metadata, datasets (IDs).
        Returns the new version.
        """
        with self.db_session_manager as session:
            current_version_id = session.query(CollectionRow.version_id).filter_by(id=collection_id).one()
            current_version = session.query(CollectionVersionRow).filter_by(version_id=current_version_id).one()
            new_version = CollectionVersionRow(**current_version.asdict)
            new_version.version_id = CollectionVersionId(self._generate_id())
            session.add(new_version)
        new_version.canonical_collection = self.get_canonical_collection(collection_id)
        return new_version

    def delete_collection_version(self, version_id: CollectionVersionId) -> None:
        """
        Deletes a collection version, if it is unpublished.
        """
        with self.db_session_manager as session:
            version = session.query(CollectionVersionRow).filter_by(version_id=version_id).one_or_none()
            if version and version.published_at is None:
                session.delete(version)

    def finalize_collection_version(
        self, collection_id: CollectionId, version_id: CollectionVersionId, published_at: Optional[datetime]
    ) -> None:
        """
        Finalizes a collection version. This is equivalent to calling:
        1. update_collection_version_mapping
        2. set_collection_version_published_at
        3. finalize_dataset_versions
        """
        self.update_collection_version_mapping(collection_id, version_id, published_at)
        self.set_collection_version_published_at(version_id, published_at)
        self.finalize_dataset_versions(version_id, published_at)

    def update_collection_version_mapping(self, collection_id: CollectionId, version_id: CollectionVersionId,
                                          published_at: datetime) -> None:
        """
        Updates the mapping between the canonical collection `collection_id` and its `version_id`
        """
        with self.db_session_manager as session:
            collection = session.query(CollectionRow).filter_by(id=collection_id).one()
            collection.version_id = version_id
            if collection.originally_published_at is None:
                collection.originally_published_at = published_at

    def set_collection_version_published_at(self, version_id: CollectionVersionId, published_at: datetime) -> None:
        """
        Sets the `published_at` datetime for a collection version
        """
        with self.db_session_manager as session:
            collection_version = session.query(CollectionVersionRow).filter_by(version_id=version_id).one()
            collection_version.published_at = published_at

    def finalize_dataset_versions(self, version_id: CollectionVersionId, published_at: datetime) -> None:
        """
        1. Updates the mapping between the canonical dataset and the latest published dataset version for each dataset
        in a CollectionVersion
        2. Sets the each canonical dataset's 'published_at' if not previously set
        """
        with self.db_session_manager as session:
            dataset_version_ids = session.query(CollectionVersionRow.datasets).filter_by(version_id=version_id).one()
            for dataset_version, dataset in (
                session.query(DatasetVersionRow, DatasetRow)
                .filter(DatasetVersionRow.dataset_id == DatasetRow.dataset_id)
                .filter(DatasetVersionRow.version_id in dataset_version_ids)
                .all()
            ):
                dataset.version_id = dataset_version.version_id
                if dataset.published_at is None:
                    dataset.published_at = published_at

    def get_dataset_version(self, dataset_version_id: DatasetVersionId) -> DatasetVersion:
        """
        Returns a dataset version by id.
        """
        with self.db_session_manager as session:
            dataset_version = session.query(DatasetVersionRow).filter_by(version_id=dataset_version_id).one()
        return self._hydrate_dataset_version(dataset_version)

    def get_all_versions_for_dataset(self, dataset_id: DatasetId) -> List[DatasetVersion]:
        """
        Returns all dataset versions for a canonical dataset_id
        """
        with self.db_session_manager as session:
            dataset_versions = session.query(DatasetVersionRow).filter_by(dataset_id=dataset_id).all()
        for i in range(len(dataset_versions)):
            dataset_versions[i] = self._hydrate_dataset_version(dataset_versions[i])
        return dataset_versions

    def get_all_datasets(self) -> Iterable[DatasetVersion]:
        """
        Returns all dataset versions.
        # TODO: Add filtering (tombstoned? remove orphaned datasets? canonical only? published?)
        """
        pass

    def get_dataset_artifacts(self, dataset_artifact_id_list: List[DatasetArtifactId]) -> List[DatasetArtifact]:
        """
        Returns all the artifacts given a list of DatasetArtifactIds
        """
        with self.db_session_manager as session:
            artifacts = session.query(DatasetArtifactRow).filter(DatasetArtifactRow.id in dataset_artifact_id_list).all() # noqa
        return artifacts

    def get_dataset_artifacts_by_version_id(self, dataset_version_id: DatasetVersionId) -> List[DatasetArtifact]:
        """
        Returns all the artifacts for a specific dataset version
        """
        with self.db_session_manager as session:
            artifact_ids = session.query(DatasetVersionRow.artifacts).filter_by(version_id=dataset_version_id).one()
        return self.get_dataset_artifacts(artifact_ids)

    def create_canonical_dataset(self, collection_version_id: CollectionVersionId) -> DatasetVersion:
        """
        Initializes a canonical dataset, generating a dataset_id and a dataset_version_id.
        Returns the newly created DatasetVersion.
        """
        with self.db_session_manager as session:
            collection_id = session.query(CollectionVersionRow.collection_id).filter_by(version_id=collection_version_id).one()
        dataset_id = DatasetId(self._generate_id())
        dataset_version_id = DatasetVersionId(self._generate_id())
        canonical_dataset = DatasetRow(dataset_id=dataset_id,
                                       dataset_version_id=dataset_version_id,
                                       published_at=None
                                       )
        dataset_version = DatasetVersionRow(version_id=dataset_version_id,
                                            dataset_id=dataset_id,
                                            collection_id=collection_id,
                                            metadata=None,
                                            artifacts=list(),
                                            status=DatasetStatus.empty(),
                                            )
        with self.db_session_manager as session:
            session.add(canonical_dataset)
            session.add(dataset_version)
        return dataset_version

    def add_dataset_artifact(self, version_id: DatasetVersionId, artifact_type: str, artifact_uri: str) -> DatasetArtifactId:
        """
        Adds a dataset artifact to an existing dataset version.
        """
        artifact_id = DatasetArtifactId(self._generate_id())
        artifact = DatasetArtifactRow(id=artifact_id,
                                      type=artifact_type,
                                      uri=artifact_uri)
        with self.db_session_manager as session:
            session.add(artifact)
            dataset_version = session.query(DatasetVersionRow).filter_by(version_id=version_id).one()
            dataset_version.artifacts.append(artifact_id)
        return artifact_id

    def update_dataset_processing_status(self, version_id: DatasetVersionId, status: DatasetProcessingStatus) -> None:
        """
        Updates the processing status for a dataset version.
        """
        with self.db_session_manager as session:
            dataset_version = session.query(DatasetVersionRow).filter_by(version_id=version_id).one()
            dataset_version.status.processing_status = status

    def update_dataset_validation_status(self, version_id: DatasetVersionId, status: DatasetValidationStatus) -> None:
        """
        Updates the validation status for a dataset version.
        """
        with self.db_session_manager as session:
            dataset_version = session.query(DatasetVersionRow).filter_by(version_id=version_id).one()
            dataset_version.status.validation_status = status

    def update_dataset_upload_status(self, version_id: DatasetVersionId, status: DatasetUploadStatus) -> None:
        """
        Updates the upload status for a dataset version.
        """
        with self.db_session_manager as session:
            dataset_version = session.query(DatasetVersionRow).filter_by(version_id=version_id).one()
            dataset_version.status.upload_status = status

    def update_dataset_conversion_status(self, version_id: DatasetVersionId, status_type: str,
                                         status: DatasetConversionStatus) -> None:
        """
        Updates the conversion status for a dataset version and for `status_type`
        """
        with self.db_session_manager as session:
            dataset_version = session.query(DatasetVersionRow).filter_by(version_id=version_id).one()
            setattr(dataset_version.status, status_type, status)

    def get_dataset_version_status(self, version_id: DatasetVersionId) -> DatasetStatus:
        """
        Returns the status for a dataset version
        """
        with self.db_session_manager as session:
            status = session.query(DatasetVersionRow.status).filter_by(version_id=version_id).one()
        return status

    def set_dataset_metadata(self, version_id: DatasetVersionId, metadata: DatasetMetadata) -> None:
        """
        Sets the metadata for a dataset version
        """
        with self.db_session_manager as session:
            dataset_version = session.query(DatasetVersionRow).filter_by(version_id=version_id).one()
            dataset_version.metadata = metadata

    def add_dataset_to_collection_version_mapping(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId
    ) -> None:
        """
        Adds a mapping between an existing collection version and a dataset version
        """
        with self.db_session_manager as session:
            collection_version = session.query(CollectionVersionRow).filter_by(version_id=collection_version_id).one()
            collection_version.datasets.append(dataset_version_id)

    def delete_dataset_from_collection_version(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId
    ) -> None:
        """
        Removes a mapping between a collection version and a dataset version
        """
        with self.db_session_manager as session:
            collection_version = session.query(CollectionVersionRow).filter_by(version_id=collection_version_id).one()
            collection_version.datasets.remove(dataset_version_id)

    def replace_dataset_in_collection_version(
        self,
        collection_version_id: CollectionVersionId,
        old_dataset_version_id: DatasetVersionId
    ) -> DatasetVersion:
        """
        Replaces an existing mapping between a collection version and a dataset version
        """
        with self.db_session_manager as session:
            collection_id = session.query(CollectionVersionRow.collection_id).filter_by(version_id=collection_version_id).one() # noqa
        new_dataset_version_id = DatasetVersionId(self._generate_id())
        new_dataset_version = DatasetVersionRow(version_id=new_dataset_version_id,
                                                dataset_id=old_dataset_version_id,
                                                collection_id=collection_id,
                                                metadata=None,
                                                artifacts=list(),
                                                status=DatasetStatus.empty()
                                                )
        with self.db_session_manager as session:
            collection_version = session.query(CollectionVersionRow).filter_by(version_id=collection_version_id).one() # noqa
            collection_version.datasets.remove(old_dataset_version_id)
            collection_version.datasets.append(new_dataset_version_id)

            canonical_dataset_id = session.query(DatasetVersionRow.dataset_id).filter_by(version_id=old_dataset_version_id).one() # noqa
            canonical_dataset = session.query(DatasetRow).filter_by(dataset_id=canonical_dataset_id).one()
        new_dataset_version.canonical_dataset = canonical_dataset
        return new_dataset_version

    def get_dataset_mapped_version(self, dataset_id: DatasetId) -> Optional[DatasetVersion]:
        """
        Returns the dataset version mapped to a canonical dataset_id, or None if not existing
        """
        with self.db_session_manager as session:
            canonical_dataset = session.query(DatasetRow).filter_by(dataset_id=dataset_id).one()
            if not canonical_dataset.version_id:
                return None
            dataset_version = session.query(DatasetVersionRow).filter_by(version_id=canonical_dataset.version_id).one()
        dataset_version.canonical_dataset = canonical_dataset
        return self._hydrate_dataset_version(dataset_version)
