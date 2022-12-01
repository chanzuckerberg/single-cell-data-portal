from datetime import datetime
from typing import Any, List, Optional, Iterable
import json
import uuid
from copy import deepcopy
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
from backend.layers.persistence.db_session import db_session_manager
from sqlalchemy import select

from backend.layers.persistence.persistence_interface import DatabaseProviderInterface


class DatabaseProvider(DatabaseProviderInterface):

    def __init__(self) -> None:
        self.db_session_manager = db_session_manager

    def _drop(self):
        from sqlalchemy.schema import DropSchema
        from backend.layers.persistence.orm import metadata
        from backend.layers.persistence.db_session import _db_session_maker
        _db_session_maker.engine.execute(DropSchema('persistence_schema', cascade=True))

    def _create(self):
        from sqlalchemy.schema import CreateSchema
        from backend.layers.persistence.orm import metadata
        from backend.layers.persistence.db_session import _db_session_maker
        _db_session_maker.engine.execute(CreateSchema('persistence_schema'))
        metadata.create_all(bind=_db_session_maker.engine)

    @staticmethod
    def _generate_id():
        return str(uuid.uuid4())

    def _row_to_collection_version(self, row: Any, canonical_collection: CanonicalCollection) -> CollectionVersion:
        return CollectionVersion(
            # self.parse_id(row.collection_id, CollectionId),
            CollectionId(str(row.collection_id)),
            CollectionVersionId(str(row.version_id)),
            row.owner,
            CollectionMetadata.from_json(row.metadata),
            row.publisher_metadata,
            [DatasetVersionId(str(dataset)) for dataset in row.datasets],
            row.published_at,
            row.created_at,
            canonical_collection
        )

    def _row_to_dataset_artifact(self, row: Any):
        pass

    def _row_to_dataset_version(self, row: Any, canonical_dataset: CanonicalDataset):
        if type(row.status) is DatasetStatus or row.status is None:
            status = row.status
        else:
            status = DatasetStatus.from_json(row.status)
        if type(row.metadata) is DatasetMetadata or row.metadata is None:
            metadata = row.metadata
        else:
            metadata = DatasetMetadata.from_json(row.metadata)
        return DatasetVersion(
            DatasetId(str(row.dataset_id)),
            DatasetVersionId(str(row.version_id)),
            CollectionId(str(row.collection_id)),
            status,
            metadata,
            row.artifacts,
            row.created_at,
            canonical_dataset
        )

    def get_canonical_collection(self, collection_id: CollectionId) -> CanonicalCollection:
        with self.db_session_manager() as session:
            collection = session.query(CollectionRow).filter_by(id=collection_id.id).one()
            return CanonicalCollection(
                CollectionId(str(collection.id)),
                collection.version_id,
                collection.originally_published_at,
                collection.tombstoned
            )

    def get_canonical_dataset(self, dataset_id: DatasetId) -> CanonicalDataset:
        with self.db_session_manager() as session:
            dataset = session.query(DatasetRow).filter_by(dataset_id=dataset_id).one()
            return CanonicalDataset(
                dataset_id,
                dataset.dataset_version_id,
                dataset.published_at
            )

    @staticmethod
    def parse_id(id: Any, id_type: Any):
        if id is None: 
            return None
        return type(id_type)(str(id))

    def create_canonical_collection(self, owner: str, collection_metadata: CollectionMetadata) -> CollectionVersion:
        """
        Creates a new canonical collection, generating a canonical collection_id and a new version_id.
        Returns the newly created CollectionVersion
        """
        collection_id = CollectionId((self._generate_id()))
        version_id = CollectionVersionId((self._generate_id()))
        now = datetime.utcnow()
        canonical_collection = CollectionRow(id=collection_id.id,
                                             version_id=version_id.id,
                                             tombstoned=False,
                                             originally_published_at=None,
                                             )

        collection_version_row = CollectionVersionRow(collection_id=collection_id.id,
                                                      version_id=version_id.id,
                                                      owner=owner,
                                                      metadata=collection_metadata.to_json(),
                                                      publisher_metadata=None,
                                                      published_at=None,
                                                      created_at=now,
                                                      datasets=list(),
                                                      )

        with self.db_session_manager() as session:
            session.add(canonical_collection)
            session.add(collection_version_row)
            return self._row_to_collection_version(collection_version_row,
                                                   CanonicalCollection(collection_id, version_id, None, False)
                                                   )

    def get_collection_version(self, version_id: CollectionVersionId) -> CollectionVersion:
        """
        Retrieves a specific collection version by id
        """
        with self.db_session_manager() as session:
            collection_version = session.query(CollectionVersionRow).filter_by(version_id=version_id.id).one_or_none()
            if collection_version is None:
                return None
            collection_id = CollectionId(str(collection_version.collection_id))
            canonical_collection = self.get_canonical_collection(collection_id)
            return self._row_to_collection_version(collection_version, canonical_collection)

    def get_collection_mapped_version(self, collection_id: CollectionId) -> Optional[CollectionVersion]:
        """
        Retrieves the latest mapped version for a collection
        """
        with self.db_session_manager() as session:
            version_id = session.query(CollectionRow.version_id).filter_by(id=collection_id.id).one_or_none()
            if version_id is None:
                return None
            collection_version = session.query(CollectionVersionRow).filter_by(version_id=version_id[0]).one()
            canonical_collection = self.get_canonical_collection(collection_id)
            return self._row_to_collection_version(collection_version, canonical_collection)

    def get_all_versions_for_collection(self, collection_id: CollectionId) -> List[CollectionVersion]:
        """
        Retrieves all versions for a specific collections, without filtering
        """
        with self.db_session_manager() as session:
            version_rows = session.query(CollectionVersionRow).filter_by(collection_id=collection_id.id).all()
            canonical_collection = self.get_canonical_collection(collection_id)
            versions = list()
            for i in range(len(version_rows)):
                version = self._row_to_collection_version(version_rows[i], canonical_collection)
                versions.append(version)
            return versions

    def get_all_collections_versions(self) -> Iterable[CollectionVersion]:
        """
        Retrieves all versions of all collections.
        TODO: for performance reasons, it might be necessary to add a filtering parameter here.
        """
        with self.db_session_manager() as session:
            versions = session.query(CollectionVersionRow).all()
            # TODO: do we need to hydrate versions with canonical collections? would require a join or many lookup calls
            return [self._row_to_collection_version(v, None) for v in versions]

    def get_all_mapped_collection_versions(self, get_tombstoned: bool = False) -> Iterable[CollectionVersion]:
        """
        Retrieves all the collection versions that are mapped to a canonical collection.
        """
        with self.db_session_manager() as session:
            if get_tombstoned:
                mapped_version_ids = session.query(CollectionRow.version_id).filter(CollectionRow.version_id.isnot(None)).all() # noqa
            else:
                mapped_version_ids = session.query(CollectionRow.version_id)\
                    .filter(CollectionRow.version_id.isnot(None))\
                    .filter_by(tombstoned=False)\
                    .all()

                # TODO: Very hacky
                mapped_version_ids = [i[0] for i in mapped_version_ids]

            versions = session.query(CollectionVersionRow).filter(CollectionVersionRow.version_id.in_(mapped_version_ids)).all() # noqa

            # TODO: do we need to hydrate versions with canonical collections? would require a join or many lookup calls
            return [self._row_to_collection_version(v, None) for v in versions]

    def delete_canonical_collection(self, collection_id: CollectionId) -> None:
        """
        Deletes (tombstones) a canonical collection.
        """
        with self.db_session_manager() as session:
            canonical_collection = session.query(CollectionRow).filter_by(id=collection_id).one()
            canonical_collection.tombstoned = True

    def save_collection_metadata(
        self, version_id: CollectionVersionId, collection_metadata: CollectionMetadata
    ) -> None:
        """
        Saves collection metadata for a collection version
        """
        with self.db_session_manager() as session:
            version = session.query(CollectionVersionRow).filter_by(version_id=version_id.id).one()
            version.metadata = collection_metadata.to_json()

    def save_collection_publisher_metadata(
            self, version_id: CollectionVersionId, publisher_metadata: Optional[dict]
    ) -> None:
        """
        Saves publisher metadata for a collection version. Specify None to remove it
        """
        print(version_id, publisher_metadata)
        with self.db_session_manager() as session:
            version = session.query(CollectionVersionRow).filter_by(version_id=version_id.id).one()
            version.publisher_metadata = publisher_metadata

    def add_collection_version(self, collection_id: CollectionId) -> CollectionVersion:
        """
        Adds a collection version to an existing canonical collection. The new version copies the following data from
         the previous version: owner, metadata, publisher_metadata, datasets (IDs).
        Returns the new version.
        """
        with self.db_session_manager() as session:
            current_version_id = session.query(CollectionRow.version_id).filter_by(id=collection_id.id).one()[0]
            current_version = session.query(CollectionVersionRow).filter_by(version_id=current_version_id).one()
            new_version = CollectionVersionRow(
                version_id=self._generate_id(),
                collection_id=collection_id.id,
                metadata=current_version.metadata,
                owner=current_version.owner,
                publisher_metadata=current_version.publisher_metadata,
                published_at=None,
                created_at=datetime.utcnow(),
                datasets=current_version.datasets
            )
            session.add(new_version)
            canonical_collection = self.get_canonical_collection(collection_id)
            return self._row_to_collection_version(new_version, canonical_collection)

    def delete_collection_version(self, version_id: CollectionVersionId) -> None:
        """
        Deletes a collection version, if it is unpublished.
        """
        with self.db_session_manager() as session:
            version = session.query(CollectionVersionRow).filter_by(version_id=version_id.id).one_or_none()
            if version and version.published_at is None:
                session.delete(version)

    def finalize_collection_version(
        self, collection_id: CollectionId, version_id: CollectionVersionId,
            published_at: Optional[datetime] = None
    ) -> None:
        """
        Finalizes a collection version. This is equivalent to calling:
        1. update_collection_version_mapping
        2. set_collection_version_published_at
        3. finalize_dataset_versions
        """
        if published_at is None:
            published_at = datetime.utcnow()
        self.update_collection_version_mapping(collection_id, version_id, published_at)
        self.set_collection_version_published_at(version_id, published_at)
        self.finalize_dataset_versions(version_id, published_at)

    def update_collection_version_mapping(self, collection_id: CollectionId, version_id: CollectionVersionId,
                                          published_at: datetime) -> None:
        """
        Updates the mapping between the canonical collection `collection_id` and its `version_id`
        """
        with self.db_session_manager() as session:
            collection = session.query(CollectionRow).filter_by(id=collection_id.id).one()
            collection.version_id = version_id.id
            if collection.originally_published_at is None:
                collection.originally_published_at = published_at

    def set_collection_version_published_at(self, version_id: CollectionVersionId, published_at: datetime) -> None:
        """
        Sets the `published_at` datetime for a collection version
        """
        with self.db_session_manager() as session:
            collection_version = session.query(CollectionVersionRow).filter_by(version_id=version_id.id).one()
            collection_version.published_at = published_at

    def finalize_dataset_versions(self, version_id: CollectionVersionId, published_at: datetime) -> None:
        """
        1. Updates the mapping between the canonical dataset and the latest published dataset version for each dataset
        in a CollectionVersion
        2. Sets the each canonical dataset's 'published_at' if not previously set
        """
        with self.db_session_manager() as session:
            dataset_version_ids = session.query(CollectionVersionRow.datasets).filter_by(version_id=version_id.id).one()[0]
            for dataset_version, dataset in (
                session.query(DatasetVersionRow, DatasetRow)
                .filter(DatasetVersionRow.dataset_id == DatasetRow.dataset_id)
                .filter(DatasetVersionRow.version_id.in_(dataset_version_ids))
                .all()
            ):
                dataset.version_id = dataset_version.version_id
                if dataset.published_at is None:
                    dataset.published_at = published_at

    def get_dataset_version(self, dataset_version_id: DatasetVersionId) -> DatasetVersion:
        """
        Returns a dataset version by id.
        """
        with self.db_session_manager() as session:
            dataset_version = session.query(DatasetVersionRow).filter_by(version_id=dataset_version_id.id).one()
            return self._row_to_dataset_version(dataset_version, self.get_canonical_dataset(dataset_version.dataset_id))

    def get_all_versions_for_dataset(self, dataset_id: DatasetId) -> List[DatasetVersion]:
        """
        Returns all dataset versions for a canonical dataset_id
        """
        dataset = self.get_canonical_dataset(dataset_id)
        with self.db_session_manager() as session:
            dataset_versions = session.query(DatasetVersionRow).filter_by(dataset_id=dataset_id.id).all()
            for i in range(len(dataset_versions)):
                dataset_versions[i] = self._row_to_dataset_version(dataset_versions[i], dataset)
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
        with self.db_session_manager() as session:
            artifacts = session.query(DatasetArtifactRow).filter(DatasetArtifactRow.id.in_(dataset_artifact_id_list)).all() # noqa
        return artifacts

    def get_dataset_artifacts_by_version_id(self, dataset_version_id: DatasetVersionId) -> List[DatasetArtifact]:
        """
        Returns all the artifacts for a specific dataset version
        """
        with self.db_session_manager() as session:
            artifact_ids = session.query(DatasetVersionRow.artifacts).filter_by(version_id=dataset_version_id).one()
        return self.get_dataset_artifacts(artifact_ids)

    def create_canonical_dataset(self, collection_version_id: CollectionVersionId) -> DatasetVersion:
        """
        Initializes a canonical dataset, generating a dataset_id and a dataset_version_id.
        Returns the newly created DatasetVersion.
        """
        with self.db_session_manager() as session:
            collection_id = session.query(CollectionVersionRow.collection_id).filter_by(version_id=collection_version_id.id).one()[0]
        dataset_id = DatasetId(self._generate_id())
        dataset_version_id = DatasetVersionId(self._generate_id())
        canonical_dataset = DatasetRow(dataset_id=dataset_id.id,
                                       dataset_version_id=dataset_version_id.id,
                                       published_at=None
                                       )
        dataset_version = DatasetVersionRow(version_id=dataset_version_id.id,
                                            dataset_id=dataset_id.id,
                                            collection_id=collection_id,
                                            metadata=None,
                                            artifacts=list(),
                                            status=DatasetStatus.empty().to_json(),
                                            created_at=datetime.utcnow(),
                                            )
        with self.db_session_manager() as session:
            session.add(canonical_dataset)
            session.add(dataset_version)
            return self._row_to_dataset_version(dataset_version, CanonicalDataset(dataset_id, dataset_version_id, None))

    def add_dataset_artifact(self, version_id: DatasetVersionId, artifact_type: str, artifact_uri: str) -> DatasetArtifactId:
        """
        Adds a dataset artifact to an existing dataset version.
        """
        artifact_id = DatasetArtifactId(self._generate_id())
        artifact = DatasetArtifactRow(id=artifact_id.id,
                                      type=artifact_type,
                                      uri=artifact_uri)
        with self.db_session_manager() as session:
            session.add(artifact)
            dataset_version = session.query(DatasetVersionRow).filter_by(version_id=version_id.id).one()
            dataset_version.artifacts.append(artifact_id.id)
        return artifact_id

    def update_dataset_processing_status(self, version_id: DatasetVersionId, status: DatasetProcessingStatus) -> None:
        """
        Updates the processing status for a dataset version.
        """
        with self.db_session_manager() as session:
            dataset_version = session.query(DatasetVersionRow).filter_by(version_id=version_id.id).one()
            dataset_version_status = json.loads(dataset_version.status)
            dataset_version_status['processing_status'] = status.value
            dataset_version.status = json.dumps(dataset_version_status)

    def update_dataset_validation_status(self, version_id: DatasetVersionId, status: DatasetValidationStatus) -> None:
        """
        Updates the validation status for a dataset version.
        """
        with self.db_session_manager() as session:
            dataset_version = session.query(DatasetVersionRow).filter_by(version_id=version_id.id).one()
            dataset_version_status = json.loads(dataset_version.status)
            dataset_version_status['validation_status'] = status.value
            dataset_version.status = json.dumps(dataset_version_status)

    def update_dataset_upload_status(self, version_id: DatasetVersionId, status: DatasetUploadStatus) -> None:
        """
        Updates the upload status for a dataset version.
        """
        with self.db_session_manager() as session:
            dataset_version = session.query(DatasetVersionRow).filter_by(version_id=version_id.id).one()
            dataset_version_status = json.loads(dataset_version.status)
            dataset_version_status['upload_status'] = status.value
            dataset_version.status = json.dumps(dataset_version_status)

    def update_dataset_conversion_status(self, version_id: DatasetVersionId, status_type: str,
                                         status: DatasetConversionStatus) -> None:
        """
        Updates the conversion status for a dataset version and for `status_type`
        """
        with self.db_session_manager() as session:
            dataset_version = session.query(DatasetVersionRow).filter_by(version_id=version_id.id).one()
            dataset_version_status = json.loads(dataset_version.status)
            dataset_version_status[status_type] = status.value
            dataset_version.status = json.dumps(dataset_version_status)

    def get_dataset_version_status(self, version_id: DatasetVersionId) -> DatasetStatus:
        """
        Returns the status for a dataset version
        """
        with self.db_session_manager() as session:
            status = session.query(DatasetVersionRow.status).filter_by(version_id=version_id.id).one()
        return json.loads(status)

    def set_dataset_metadata(self, version_id: DatasetVersionId, metadata: DatasetMetadata) -> None:
        """
        Sets the metadata for a dataset version
        """
        with self.db_session_manager() as session:
            dataset_version = session.query(DatasetVersionRow).filter_by(version_id=version_id.id).one()
            dataset_version.metadata = metadata.to_json()

    def add_dataset_to_collection_version_mapping(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId
    ) -> None:
        """
        Adds a mapping between an existing collection version and a dataset version
        """
        with self.db_session_manager() as session:
            collection_version = session.query(CollectionVersionRow).filter_by(version_id=collection_version_id.id).one()
            datasets = [uuid.UUID(dataset_version_id.id)]
            datasets.extend(collection_version.datasets)
            collection_version.datasets = datasets

    def delete_dataset_from_collection_version(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId
    ) -> None:
        """
        Removes a mapping between a collection version and a dataset version
        """
        with self.db_session_manager() as session:
            collection_version = session.query(CollectionVersionRow).filter_by(version_id=collection_version_id.id).one()
            datasets = deepcopy(collection_version.datasets)
            datasets.remove(uuid.UUID(dataset_version_id.id))
            collection_version.datasets = datasets

    def replace_dataset_in_collection_version(
        self,
        collection_version_id: CollectionVersionId,
        old_dataset_version_id: DatasetVersionId
    ) -> DatasetVersion:
        """
        Replaces an existing mapping between a collection version and a dataset version
        """
        with self.db_session_manager() as session:
            collection_id = session.query(CollectionVersionRow.collection_id).filter_by(version_id=collection_version_id.id).one()[0] # noqa
            dataset_id = session.query(DatasetVersionRow.dataset_id).filter_by(version_id=old_dataset_version_id.id).one()[0]
            new_dataset_version_id = DatasetVersionId(self._generate_id())
            new_dataset_version = DatasetVersionRow(version_id=new_dataset_version_id.id,
                                                    dataset_id=dataset_id,
                                                    collection_id=collection_id,
                                                    metadata=None,
                                                    artifacts=list(),
                                                    status=DatasetStatus.empty().to_json(),
                                                    created_at=datetime.utcnow(),
                                                    )
            session.add(new_dataset_version)
            canonical_dataset = self.get_canonical_dataset(dataset_id)
            self.add_dataset_to_collection_version_mapping(collection_version_id, new_dataset_version_id)
            self.delete_dataset_from_collection_version(collection_version_id, old_dataset_version_id)
            return self._row_to_dataset_version(new_dataset_version, canonical_dataset)

    def get_dataset_mapped_version(self, dataset_id: DatasetId) -> Optional[DatasetVersion]:
        """
        Returns the dataset version mapped to a canonical dataset_id, or None if not existing
        """
        with self.db_session_manager() as session:
            canonical_dataset = session.query(DatasetRow).filter_by(dataset_id=dataset_id).one()
            if not canonical_dataset.version_id:
                return None
            dataset_version = session.query(DatasetVersionRow).filter_by(version_id=canonical_dataset.version_id).one()
        dataset_version.canonical_dataset = canonical_dataset
        return self._hydrate_dataset_version(dataset_version)
