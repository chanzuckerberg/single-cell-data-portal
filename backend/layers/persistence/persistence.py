import contextlib
import json
import logging
import uuid
from contextlib import contextmanager
from datetime import datetime
from typing import Any, Iterable, List, Optional, Tuple, Union

from sqlalchemy import create_engine, delete
from sqlalchemy.exc import ProgrammingError, SQLAlchemyError
from sqlalchemy.orm import Session, sessionmaker

from backend.common.corpora_config import CorporaDbConfig
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
    DatasetArtifactType,
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
from backend.layers.common.helpers import set_revised_at_field
from backend.layers.persistence.constants import SCHEMA_NAME
from backend.layers.persistence.orm import (
    CollectionTable,
    CollectionVersionTable,
    DatasetArtifactTable,
    DatasetTable,
    DatasetVersionTable,
)
from backend.layers.persistence.persistence_interface import DatabaseProviderInterface, PersistenceException

logger = logging.getLogger(__name__)


class DatabaseProvider(DatabaseProviderInterface):
    def __init__(self, database_uri: str = None, schema_name: str = SCHEMA_NAME) -> None:
        if not database_uri:
            database_uri = CorporaDbConfig().database_uri
        self._engine = create_engine(database_uri, connect_args={"connect_timeout": 5})
        self._session_maker = sessionmaker(bind=self._engine)
        self._schema_name = schema_name
        with contextlib.suppress(Exception):
            self._create_schema()

    def _drop_schema(self):
        from sqlalchemy.schema import DropSchema

        with contextlib.suppress(ProgrammingError):
            self._engine.execute(DropSchema(self._schema_name, cascade=True))

    def _create_schema(self):
        from sqlalchemy.schema import CreateSchema

        from backend.layers.persistence.orm import metadata

        self._engine.execute(CreateSchema(self._schema_name))
        metadata.schema = self._schema_name
        metadata.create_all(bind=self._engine)

    @contextmanager
    def _manage_session(self, **kwargs):
        try:
            if not self._session_maker:
                self._session_maker = sessionmaker(bind=self._engine)
            session = self._session_maker(**kwargs)
            yield session
            if session.transaction:
                session.commit()
            else:
                session.expire_all()
        except SQLAlchemyError as e:
            logger.exception(e)
            if session is not None:
                session.rollback()
            raise PersistenceException("Failed to commit.") from None
        finally:
            if session is not None:
                session.close()

    def _row_to_collection_version(self, row: Any, canonical_collection: CanonicalCollection) -> CollectionVersion:
        return CollectionVersion(
            collection_id=CollectionId(str(row.collection_id)),
            version_id=CollectionVersionId(str(row.id)),
            owner=row.owner,
            curator_name=row.curator_name,
            metadata=CollectionMetadata.from_json(row.collection_metadata),
            publisher_metadata=None if row.publisher_metadata is None else json.loads(row.publisher_metadata),
            datasets=[DatasetVersionId(str(id)) for id in row.datasets],
            published_at=row.published_at,
            created_at=row.created_at,
            schema_version=row.schema_version,
            canonical_collection=canonical_collection,
        )

    def _row_to_collection_version_with_datasets(
        self, row: Any, canonical_collection: CanonicalCollection, datasets: List[DatasetVersion]
    ) -> CollectionVersionWithDatasets:
        return CollectionVersionWithDatasets(
            collection_id=CollectionId(str(row.collection_id)),
            version_id=CollectionVersionId(str(row.id)),
            owner=row.owner,
            curator_name=row.curator_name,
            metadata=CollectionMetadata.from_json(row.collection_metadata),
            publisher_metadata=None if row.publisher_metadata is None else json.loads(row.publisher_metadata),
            datasets=datasets,
            published_at=row.published_at,
            created_at=row.created_at,
            schema_version=row.schema_version,
            canonical_collection=canonical_collection,
        )

    def _row_to_canonical_dataset(self, row: Any):
        return CanonicalDataset(
            DatasetId(str(row.id)),
            None if row.version_id is None else DatasetVersionId(str(row.version_id)),
            row.tombstone,
            row.published_at,
        )

    def _row_to_dataset_artifact(self, row: Any):
        return DatasetArtifact(
            DatasetArtifactId(str(row.id)),
            row.type,
            row.uri,
        )

    def _row_to_dataset_version(self, row: Any, canonical_dataset: CanonicalDataset, artifacts: List[DatasetArtifact]):
        if type(row.status) is DatasetStatus or row.status is None:
            status = row.status
        else:
            status = DatasetStatus.from_json(row.status)
        if type(row.dataset_metadata) is DatasetMetadata or row.dataset_metadata is None:
            metadata = row.dataset_metadata
        else:
            metadata = DatasetMetadata.from_json(row.dataset_metadata)
        return DatasetVersion(
            DatasetId(str(row.dataset_id)),
            DatasetVersionId(str(row.id)),
            CollectionId(str(row.collection_id)),
            status,
            metadata,
            artifacts,
            row.created_at,
            canonical_dataset,
        )

    def _hydrate_dataset_version(self, dataset_version: DatasetVersionTable) -> DatasetVersion:
        """
        Populates canonical_dataset, artifacts, and status for single DatasetVersionRow
        """
        canonical_dataset = self.get_canonical_dataset(DatasetId(str(dataset_version.dataset_id)))
        artifacts = self.get_dataset_artifacts(dataset_version.artifacts)
        return self._row_to_dataset_version(dataset_version, canonical_dataset, artifacts)

    def get_canonical_collection(self, collection_id: CollectionId) -> CanonicalCollection:
        with self._manage_session() as session:
            collection = session.query(CollectionTable).filter_by(id=collection_id.id).one_or_none()
            if collection is None:
                return None
            return CanonicalCollection(
                CollectionId(str(collection.id)),
                None if collection.version_id is None else CollectionVersionId(str(collection.version_id)),
                collection.originally_published_at,
                collection.revised_at,
                collection.tombstone,
            )

    def get_canonical_dataset(self, dataset_id: DatasetId) -> CanonicalDataset:
        with self._manage_session() as session:
            dataset = session.query(DatasetTable).filter_by(id=dataset_id.id).one_or_none()
            if dataset is None:
                return None
            return CanonicalDataset(
                dataset_id,
                None if dataset.version_id is None else DatasetVersionId(str(dataset.version_id)),
                dataset.tombstone,
                dataset.published_at,
            )

    def create_canonical_collection(
        self, owner: str, curator_name: str, collection_metadata: CollectionMetadata
    ) -> CollectionVersion:
        """
        Creates a new canonical collection, generating a canonical collection_id and a new version_id.
        Returns the newly created CollectionVersion
        """
        collection_id = CollectionId()
        version_id = CollectionVersionId()
        now = datetime.utcnow()
        canonical_collection = CollectionTable(
            id=collection_id.id, version_id=None, tombstone=False, originally_published_at=None, revised_at=None
        )

        collection_version_row = CollectionVersionTable(
            id=version_id.id,
            collection_id=collection_id.id,
            owner=owner,
            curator_name=curator_name,
            collection_metadata=collection_metadata.to_json(),
            publisher_metadata=None,
            published_at=None,
            created_at=now,
            schema_version=None,
            datasets=list(),
        )

        with self._manage_session() as session:
            session.add(canonical_collection)
            session.add(collection_version_row)

            return self._row_to_collection_version(
                collection_version_row, CanonicalCollection(collection_id, None, None, None, False)
            )

    def get_collection_version(self, version_id: CollectionVersionId) -> CollectionVersion:
        """
        Retrieves a specific collection version by id
        """
        with self._manage_session() as session:
            collection_version = session.query(CollectionVersionTable).filter_by(id=version_id.id).one_or_none()
            if collection_version is None:
                return None
            collection_id = CollectionId(str(collection_version.collection_id))
            canonical_collection = self.get_canonical_collection(collection_id)
            return self._row_to_collection_version(collection_version, canonical_collection)

    def get_dataset_versions_by_id(
        self, ids: List[DatasetVersionId], get_tombstoned: bool = False
    ) -> List[DatasetVersion]:
        ids = [dv_id.id for dv_id in ids]
        with self._manage_session() as session:
            versions = session.query(DatasetVersionTable).filter(DatasetVersionTable.id.in_(ids)).all()
            canonical_ids = []
            artifact_ids = []
            for version in versions:
                canonical_ids.append(version.dataset_id)
                artifact_ids.extend(version.artifacts)

            if get_tombstoned:
                canonical_dataset_query = session.query(DatasetTable).filter(DatasetTable.id.in_(canonical_ids))
            else:
                canonical_dataset_query = (
                    session.query(DatasetTable)
                    .filter(DatasetTable.id.in_(canonical_ids))
                    .filter(DatasetTable.tombstone.is_(False))
                )
            canonical_datasets = canonical_dataset_query.all()

            canonical_map = {canonical_dataset.id: canonical_dataset for canonical_dataset in canonical_datasets}

            artifacts = session.query(DatasetArtifactTable).filter(DatasetArtifactTable.id.in_(artifact_ids)).all()
            artifact_map = {artifact.id: artifact for artifact in artifacts}

            datasets = []
            for version in versions:
                canonical_dataset_row = canonical_map.get(version.dataset_id)
                if not canonical_dataset_row:
                    continue  # Dataset has the wrong tombstone value
                canonical_dataset = self._row_to_canonical_dataset(canonical_dataset_row)
                version_artifacts = [
                    self._row_to_dataset_artifact(artifact_map.get(artifact_id)) for artifact_id in version.artifacts
                ]
                datasets.append(self._row_to_dataset_version(version, canonical_dataset, version_artifacts))
        return datasets

    def get_collection_version_with_datasets(
        self, version_id: CollectionVersionId, get_tombstoned: bool = False
    ) -> CollectionVersionWithDatasets:
        """
        Retrieves a specific collection version by id, with datasets
        """
        with self._manage_session() as session:
            collection_version = session.query(CollectionVersionTable).filter_by(id=version_id.id).one_or_none()
            if collection_version is None:
                return None
            if not get_tombstoned:
                collection_exists = (
                    session.query(CollectionTable.id)
                    .filter_by(id=collection_version.collection_id, tombstone=False)
                    .one_or_none()
                )
                if not collection_exists:
                    return None
            collection_id = CollectionId(str(collection_version.collection_id))
            canonical_collection = self.get_canonical_collection(collection_id)
            all_collection_versions_rows = (
                session.query(CollectionVersionTable).filter_by(collection_id=canonical_collection.id.id).all()
            )
            all_collection_versions = [
                self._row_to_collection_version(c_v_row, canonical_collection)
                for c_v_row in all_collection_versions_rows
            ]
            dataset_versions = self.get_dataset_versions_by_id(
                [DatasetVersionId(str(id)) for id in collection_version.datasets]
            )
            set_revised_at_field(dataset_versions, all_collection_versions)
            return self._row_to_collection_version_with_datasets(
                collection_version, canonical_collection, dataset_versions
            )

    def get_collection_mapped_version(self, collection_id: CollectionId) -> Optional[CollectionVersionWithDatasets]:
        """
        Retrieves the latest mapped version for a collection
        """
        with self._manage_session() as session:
            version_id = session.query(CollectionTable.version_id).filter_by(id=collection_id.id).one_or_none()
            # TODO: figure out this hack
            if version_id is None or version_id[0] is None:
                return None
            version_id = version_id[0]
            collection_versions = session.query(CollectionVersionTable).filter_by(collection_id=collection_id.id).all()

            collection_version = next(c_v_row for c_v_row in collection_versions if c_v_row.id == version_id)
            canonical_collection = self.get_canonical_collection(collection_id)
            dataset_versions = self.get_dataset_versions_by_id(
                [DatasetVersionId(str(id)) for id in collection_version.datasets]
            )
            all_collection_versions = [
                self._row_to_collection_version(c_v_row, canonical_collection) for c_v_row in collection_versions
            ]
            set_revised_at_field(dataset_versions, all_collection_versions)
            return self._row_to_collection_version_with_datasets(
                collection_version, canonical_collection, dataset_versions
            )

    def get_all_versions_for_collection(self, collection_id: CollectionId) -> List[CollectionVersionWithDatasets]:
        """
        Retrieves all versions for a specific collections, without filtering
        """
        with self._manage_session() as session:
            version_rows = session.query(CollectionVersionTable).filter_by(collection_id=collection_id.id).all()
            canonical_collection = self.get_canonical_collection(collection_id)
            versions = list()
            for i in range(len(version_rows)):
                datasets = self.get_dataset_versions_by_id(
                    [DatasetVersionId(str(id)) for id in version_rows[i].datasets]
                )
                version = self._row_to_collection_version_with_datasets(version_rows[i], canonical_collection, datasets)
                versions.append(version)
            return versions

    def get_all_collections_versions(self, get_tombstoned: bool = False) -> Iterable[CollectionVersion]:
        """
        Retrieves all versions of all collections.
        TODO: for performance reasons, it might be necessary to add a filtering parameter here.
        """
        with self._manage_session() as session:
            versions = session.query(CollectionVersionTable).all()

            # Create a canonical mapping
            if get_tombstoned:
                all_canonical_collections = session.query(CollectionTable)
            else:
                all_canonical_collections = session.query(CollectionTable).filter(CollectionTable.tombstone.isnot(True))

            all_canonical_map = dict()
            for collection_row in all_canonical_collections.all():
                all_canonical_map[str(collection_row.id)] = CanonicalCollection(
                    CollectionId(str(collection_row.id)),
                    None if collection_row.version_id is None else CollectionVersionId(str(collection_row.version_id)),
                    collection_row.originally_published_at,
                    collection_row.revised_at,
                    collection_row.tombstone,
                )

            result = []
            all_dataset_tombstones = {
                str(dataset.id)
                for dataset in session.query(DatasetTable).filter(DatasetTable.tombstone.is_(True)).all()
            }
            all_dataset_version_mappings = {
                str(dataset_version.id): str(dataset_version.dataset_id)
                for dataset_version in session.query(DatasetVersionTable).all()
            }
            for v in versions:
                include_dataset_version_ids = []
                if str(v.collection_id) in all_canonical_map:
                    for dataset_version_id in v.datasets:
                        dataset_version_id_str = str(dataset_version_id)
                        dataset_id = all_dataset_version_mappings[dataset_version_id_str]
                        if dataset_id:
                            if not get_tombstoned and dataset_id in all_dataset_tombstones:
                                continue
                            include_dataset_version_ids.append(dataset_version_id)
                    v.datasets = include_dataset_version_ids
                    result.append(self._row_to_collection_version(v, all_canonical_map[str(v.collection_id)]))

            return result

    def get_all_mapped_collection_versions(self, get_tombstoned: bool = False) -> Iterable[CollectionVersion]:
        """
        Retrieves all the collection versions that are mapped to a canonical collection. This method does not require
        tombstone filtering at the 'individual Dataset' level because, by definition, only active Dataset version ids
        will be present in the CollectionVersion.datasets array for active (mapped) Collection versions.
        """
        with self._manage_session() as session:
            if get_tombstoned:
                canonical_collections = session.query(CollectionTable).filter(CollectionTable.version_id.isnot(None))
            else:
                canonical_collections = (
                    session.query(CollectionTable)
                    .filter(CollectionTable.version_id.isnot(None))
                    .filter_by(tombstone=False)
                )

            mapped_version_ids = {cc.version_id: cc for cc in canonical_collections.all()}
            versions = (
                session.query(CollectionVersionTable)
                .filter(CollectionVersionTable.id.in_(mapped_version_ids.keys()))
                .all()
            )  # noqa

            for version in versions:
                canonical_row = mapped_version_ids[version.id]
                canonical = CanonicalCollection(
                    CollectionId(str(canonical_row.id)),
                    CollectionVersionId(str(canonical_row.version_id)),
                    canonical_row.originally_published_at,
                    canonical_row.revised_at,
                    canonical_row.tombstone,
                )
                yield self._row_to_collection_version(version, canonical)

    def tombstone_collection(self, collection_id: CollectionId) -> None:
        """
        Deletes (tombstones) a canonical collection.
        """
        with self._manage_session() as session:
            canonical_collection = session.query(CollectionTable).filter_by(id=collection_id.id).one_or_none()
            if canonical_collection:
                canonical_collection.tombstone = True
            # Tombstone the Datasets individually as well; technically not necessary but will protect us from bugs
            dataset_versions = session.query(DatasetVersionTable).filter_by(collection_id=collection_id.id).all()
            datasets = session.query(DatasetTable).filter(
                DatasetTable.id.in_([dv.dataset_id for dv in dataset_versions])
            )
            for dataset in datasets:
                dataset.tombstone = True

    def save_collection_metadata(
        self, version_id: CollectionVersionId, collection_metadata: CollectionMetadata
    ) -> None:
        """
        Saves collection metadata for a collection version
        """
        with self._manage_session() as session:
            version = session.query(CollectionVersionTable).filter_by(id=version_id.id).one()
            version.collection_metadata = collection_metadata.to_json()

    def save_collection_publisher_metadata(
        self, version_id: CollectionVersionId, publisher_metadata: Optional[dict]
    ) -> None:
        """
        Saves publisher metadata for a collection version. Specify None to remove it
        """
        with self._manage_session() as session:
            version = session.query(CollectionVersionTable).filter_by(id=version_id.id).one()
            version.publisher_metadata = json.dumps(publisher_metadata)

    def add_collection_version(self, collection_id: CollectionId) -> CollectionVersionId:
        """
        Adds a collection version to an existing canonical collection. The new version copies all data from
        the previous version except version_id, schema_version, and datetime-based fields (i.e. created_at,
        published_at)
        Returns the new version id.
        """
        with self._manage_session() as session:
            current_version_id = session.query(CollectionTable.version_id).filter_by(id=collection_id.id).one()[0]
            current_version = session.query(CollectionVersionTable).filter_by(id=current_version_id).one()
            new_version_id = CollectionVersionId()
            new_version = CollectionVersionTable(
                id=new_version_id.id,
                collection_id=collection_id.id,
                collection_metadata=current_version.collection_metadata,
                owner=current_version.owner,
                curator_name=current_version.curator_name,
                publisher_metadata=current_version.publisher_metadata,
                published_at=None,
                created_at=datetime.utcnow(),
                schema_version=None,
                datasets=current_version.datasets,
            )
            session.add(new_version)
            return CollectionVersionId(new_version_id)

    def delete_collection(self, collection_id: CollectionId) -> None:
        """
        Delete an unpublished Collection
        """
        with self._manage_session() as session:
            collection = session.query(CollectionTable).filter_by(id=collection_id.id).one_or_none()
            if collection:
                if collection.originally_published_at:
                    raise CollectionIsPublishedException(f"Published Collection {collection_id} cannot be deleted")
                session.delete(collection)

    def delete_collection_version(self, version_id: CollectionVersionId) -> None:
        """
        Deletes a collection version, if it is unpublished.
        """
        with self._manage_session() as session:
            version = session.query(CollectionVersionTable).filter_by(id=version_id.id).one_or_none()
            if version:
                if version.published_at:
                    raise CollectionIsPublishedException(f"Published Collection Version {version_id} cannot be deleted")
                session.delete(version)

    def delete_datasets(self, datasets: List[Union[DatasetId, CanonicalDataset]]) -> None:
        """
        Delete an unpublished DatasetTable row (and its dependent DatasetVersionTable and DatasetArtifactTable rows)
        """
        with self._manage_session() as session:
            for d in datasets:
                d_id = d.id if isinstance(d, DatasetId) else d.dataset_id.id
                dataset_row = session.query(DatasetTable).filter_by(id=d_id).one()
                if dataset_row.published_at:
                    raise DatasetIsPublishedException(f"Published Dataset {d_id} cannot be deleted")
                dataset_versions = session.query(DatasetVersionTable).filter_by(dataset_id=d_id).all()
                self._delete_dataset_version_and_artifact_rows(dataset_versions, session)
                session.delete(dataset_row)

    def delete_dataset_versions(self, dataset_versions: List[Union[DatasetVersionId, DatasetVersion]]) -> None:
        """
        Deletes DatasetVersionTable rows.
        """
        with self._manage_session() as session:
            ids = [
                str(d_v.id) if isinstance(d_v, DatasetVersionId) else str(d_v.version_id.id) for d_v in dataset_versions
            ]
            dataset_version_rows = session.query(DatasetVersionTable).filter(DatasetVersionTable.id.in_(ids)).all()
            self._delete_dataset_version_and_artifact_rows(dataset_version_rows, session)

    def _delete_dataset_version_and_artifact_rows(
        self, dataset_version_rows: List[DatasetVersionTable], session: Session
    ) -> None:
        """
        Delete DatasetVersionTable rows (and their dependent DatasetArtifactTable rows)
        """
        for d_v_row in dataset_version_rows:
            ids = [str(_id) for _id in d_v_row.artifacts]
            artifact_delete_statement = delete(DatasetArtifactTable).where(DatasetArtifactTable.id.in_(ids))
            session.execute(artifact_delete_statement)
            session.delete(d_v_row)
        session.flush()

    def finalize_collection_version(
        self,
        collection_id: CollectionId,
        version_id: CollectionVersionId,
        schema_version: str,
        published_at: Optional[datetime] = None,
        update_revised_at: bool = False,
    ) -> List[str]:
        """
        Finalizes a collection version. Returns a list of ids for all Dataset Versions for any/all tombstoned Datasets.
        """
        published_at = published_at if published_at else datetime.utcnow()
        with self._manage_session() as session:
            # update canonical collection -> collection version mapping
            collection = session.query(CollectionTable).filter_by(id=collection_id.id).one()
            previous_c_v_id = collection.version_id

            collection.version_id = version_id.id
            # update canonical collection timestamps depending on whether this is its first publish
            if collection.originally_published_at is None:
                collection.originally_published_at = published_at
            # if not first publish, update revised_at if flagged to do so
            elif update_revised_at:
                collection.revised_at = published_at

            # update collection version
            collection_version = session.query(CollectionVersionTable).filter_by(id=version_id.id).one()
            collection_version.published_at = published_at
            collection_version.schema_version = schema_version

            dataset_ids_for_new_collection_version = [
                d.dataset_id.id for d in self.get_collection_version_with_datasets(version_id).datasets
            ]

            # finalize collection version's dataset versions
            dataset_ids_to_tombstone = []
            if previous_c_v_id:
                # Publishing a revision; Datasets could have been removed
                previous_collection_version = self.get_collection_version_with_datasets(
                    CollectionVersionId(previous_c_v_id)
                )
                previous_d_ids = [d.dataset_id.id for d in previous_collection_version.datasets]
                for previous_d_id in previous_d_ids:
                    if previous_d_id not in dataset_ids_for_new_collection_version:
                        dataset_ids_to_tombstone.append(previous_d_id)

            # get all dataset versions for the datasets that are being tombstoned
            dataset_version_ids_to_delete_from_s3 = []
            if dataset_ids_to_tombstone:
                datasets = session.query(DatasetTable).filter(DatasetTable.id.in_(dataset_ids_to_tombstone)).all()
                for dataset in datasets:
                    dataset.tombstone = True
                    dataset_all_versions = session.query(DatasetVersionTable).filter_by(dataset_id=dataset.id).all()
                    dataset_version_ids_to_delete_from_s3.extend([dv.id for dv in dataset_all_versions])

            # update dataset versions for datasets that are not being tombstoned
            dataset_version_ids = session.query(CollectionVersionTable.datasets).filter_by(id=version_id.id).one()[0]
            for dataset_version, dataset in (
                session.query(DatasetVersionTable, DatasetTable)
                .filter(DatasetVersionTable.dataset_id == DatasetTable.id)
                .filter(DatasetVersionTable.id.in_(dataset_version_ids))
                .all()
            ):
                dataset.version_id = dataset_version.id  # point the dataset to the new version
                if dataset.published_at is None:
                    dataset.published_at = published_at

            return dataset_version_ids_to_delete_from_s3

    def get_dataset_version(self, dataset_version_id: DatasetVersionId, get_tombstoned: bool = False) -> DatasetVersion:
        """
        Returns a dataset version by id.
        """
        with self._manage_session() as session:
            dataset_version = session.query(DatasetVersionTable).filter_by(id=dataset_version_id.id).one_or_none()
            if dataset_version is None:
                return None
            if not get_tombstoned:
                dataset_exists = (
                    session.query(DatasetTable.id)
                    .filter_by(id=dataset_version.dataset_id, tombstone=False)
                    .one_or_none()
                )
                if not dataset_exists:
                    return None
            return self._hydrate_dataset_version(dataset_version)

    def get_all_dataset_versions_for_collection(
        self, collection_id: CollectionId, from_date: datetime = datetime.min
    ) -> List[DatasetVersion]:
        """
        Get all Dataset versions -- published and unpublished -- for a canonical Collection
        """
        from_date = datetime.min if from_date is None else from_date
        with self._manage_session() as session:
            dataset_versions = (
                session.query(DatasetVersionTable)
                .filter(DatasetVersionTable.collection_id == uuid.UUID(collection_id.id))
                .filter(DatasetVersionTable.created_at >= from_date)
                .all()
            )
            return [self._hydrate_dataset_version(dv) for dv in dataset_versions]

    def get_all_versions_for_dataset(self, dataset_id: DatasetId) -> List[DatasetVersion]:
        """
        Returns all dataset versions for a canonical dataset_id
        """
        dataset = self.get_canonical_dataset(dataset_id)
        with self._manage_session() as session:
            dataset_versions = session.query(DatasetVersionTable).filter_by(dataset_id=dataset_id.id).all()
            artifact_ids = [artifact_id for dv in dataset_versions for artifact_id in dv.artifacts]
            artifacts = session.query(DatasetArtifactTable).filter(DatasetArtifactTable.id.in_(artifact_ids)).all()
            artifact_map = {artifact.id: artifact for artifact in artifacts}
            for i in range(len(dataset_versions)):
                version = dataset_versions[i]
                version_artifacts = [
                    self._row_to_dataset_artifact(artifact_map.get(artifact_id)) for artifact_id in version.artifacts
                ]
                dataset_versions[i] = self._row_to_dataset_version(version, dataset, version_artifacts)
            return dataset_versions

    def get_all_mapped_datasets_and_collections(self) -> Tuple[List[DatasetVersion], List[CollectionVersion]]:
        """
        Returns all mapped datasets and mapped collection versions.
        """
        active_collections = list(self.get_all_mapped_collection_versions())
        dataset_version_ids = []
        for collection in active_collections:
            dataset_version_ids.extend(collection.datasets)
        return list(self.get_dataset_versions_by_id(dataset_version_ids)), active_collections

    def get_dataset_artifacts(self, dataset_artifact_id_list: List[DatasetArtifactId]) -> List[DatasetArtifact]:
        """
        Returns all the artifacts given a list of DatasetArtifactIds
        """
        with self._manage_session() as session:
            artifacts = (
                session.query(DatasetArtifactTable)
                .filter(DatasetArtifactTable.id.in_([str(i) for i in dataset_artifact_id_list]))
                .all()
            )  # noqa
            return [self._row_to_dataset_artifact(a) for a in artifacts]

    def get_dataset_artifacts_by_version_id(self, dataset_version_id: DatasetVersionId) -> List[DatasetArtifact]:
        """
        Returns all the artifacts for a specific dataset version
        """
        with self._manage_session() as session:
            artifact_ids = (
                session.query(DatasetVersionTable.artifacts).filter_by(version_id=dataset_version_id.id).one()
            )
        return self.get_dataset_artifacts(artifact_ids[0])

    def create_canonical_dataset(self, collection_version_id: CollectionVersionId) -> DatasetVersion:
        """
        Initializes a canonical dataset, generating a dataset_id and a dataset_version_id.
        Returns the newly created DatasetVersion.
        """
        with self._manage_session() as session:
            collection_id = (
                session.query(CollectionVersionTable.collection_id).filter_by(id=collection_version_id.id).one()[0]
            )
        dataset_id = DatasetId()
        dataset_version_id = DatasetVersionId()
        canonical_dataset = DatasetTable(id=dataset_id.id, version_id=None, published_at=None, tombstone=False)
        dataset_version = DatasetVersionTable(
            id=dataset_version_id.id,
            dataset_id=dataset_id.id,
            collection_id=collection_id,
            dataset_metadata=None,
            artifacts=list(),
            status=DatasetStatus.empty().to_json(),
            created_at=datetime.utcnow(),
        )

        with self._manage_session() as session:
            session.add(canonical_dataset)
            session.add(dataset_version)
            return self._row_to_dataset_version(dataset_version, CanonicalDataset(dataset_id, None, False, None), [])

    def add_dataset_artifact(
        self, version_id: DatasetVersionId, artifact_type: DatasetArtifactType, artifact_uri: str
    ) -> DatasetArtifactId:
        """
        Adds a dataset artifact to an existing dataset version.
        """
        artifact_id = DatasetArtifactId()
        artifact = DatasetArtifactTable(id=artifact_id.id, type=artifact_type, uri=artifact_uri)
        with self._manage_session() as session:
            session.add(artifact)
            dataset_version = session.query(DatasetVersionTable).filter_by(id=version_id.id).one()
            artifacts = list(dataset_version.artifacts)
            artifacts.append(uuid.UUID(artifact_id.id))
            dataset_version.artifacts = artifacts
        return artifact_id

    def update_dataset_artifact(self, artifact_id: DatasetArtifactId, artifact_uri: str) -> None:
        """
        Updates uri for an existing artifact_id
        """
        with self._manage_session() as session:
            artifact = session.query(DatasetArtifactTable).filter_by(id=artifact_id.id).one()
            artifact.uri = artifact_uri

    def update_dataset_processing_status(self, version_id: DatasetVersionId, status: DatasetProcessingStatus) -> None:
        """
        Updates the processing status for a dataset version.
        """
        with self._manage_session() as session:
            dataset_version = session.query(DatasetVersionTable).filter_by(id=version_id.id).one()
            dataset_version_status = json.loads(dataset_version.status)
            dataset_version_status["processing_status"] = status.value
            dataset_version.status = json.dumps(dataset_version_status)

    def update_dataset_validation_status(self, version_id: DatasetVersionId, status: DatasetValidationStatus) -> None:
        """
        Updates the validation status for a dataset version.
        """
        with self._manage_session() as session:
            dataset_version = session.query(DatasetVersionTable).filter_by(id=version_id.id).one()
            dataset_version_status = json.loads(dataset_version.status)
            dataset_version_status["validation_status"] = status.value
            dataset_version.status = json.dumps(dataset_version_status)

    def update_dataset_upload_status(self, version_id: DatasetVersionId, status: DatasetUploadStatus) -> None:
        """
        Updates the upload status for a dataset version.
        """
        with self._manage_session() as session:
            dataset_version = session.query(DatasetVersionTable).filter_by(id=version_id.id).one()
            dataset_version_status = json.loads(dataset_version.status)
            dataset_version_status["upload_status"] = status.value
            dataset_version.status = json.dumps(dataset_version_status)

    def update_dataset_conversion_status(
        self, version_id: DatasetVersionId, status_type: str, status: DatasetConversionStatus
    ) -> None:
        """
        Updates the conversion status for a dataset version and for `status_type`
        """
        with self._manage_session() as session:
            dataset_version = session.query(DatasetVersionTable).filter_by(id=version_id.id).one()
            dataset_version_status = json.loads(dataset_version.status)
            dataset_version_status[status_type] = status.value
            dataset_version.status = json.dumps(dataset_version_status)

    def update_dataset_validation_message(self, version_id: DatasetVersionId, validation_message: str) -> None:
        with self._manage_session() as session:
            dataset_version = session.query(DatasetVersionTable).filter_by(id=version_id.id).one()
            dataset_version_status = json.loads(dataset_version.status)
            dataset_version_status["validation_message"] = validation_message
            dataset_version.status = json.dumps(dataset_version_status)

    def get_dataset_version_status(self, version_id: DatasetVersionId) -> DatasetStatus:
        """
        Returns the status for a dataset version
        """
        with self._manage_session() as session:
            status = session.query(DatasetVersionTable.status).filter_by(id=version_id.id).one()
        return DatasetStatus.from_json(status[0])

    def set_dataset_metadata(self, version_id: DatasetVersionId, metadata: DatasetMetadata) -> None:
        """
        Sets the metadata for a dataset version
        """
        with self._manage_session() as session:
            dataset_version = session.query(DatasetVersionTable).filter_by(id=version_id.id).one()
            dataset_version.dataset_metadata = metadata.to_json()

    def add_dataset_to_collection_version_mapping(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId
    ) -> None:
        """
        Adds a mapping between an existing collection version and a dataset version
        """
        with self._manage_session() as session:
            collection_version = session.query(CollectionVersionTable).filter_by(id=collection_version_id.id).one()
            # TODO: alternatively use postgres `array_append`
            # TODO: make sure that the UUID conversion works
            updated_datasets = list(collection_version.datasets)
            updated_datasets.append(uuid.UUID(dataset_version_id.id))
            collection_version.datasets = updated_datasets

    def delete_dataset_from_collection_version(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId
    ) -> None:
        """
        Removes a mapping between a collection version and a dataset version.

        :param collection_version_id: the CollectionVersion or the CollectionVersionId
        :param dataset_version_id: the DatasetVersionId
        :param delete_dv_row: boolean flag - when True, delete DatasetVersion row (and dependent DatasetArtifact rows)
        """
        with self._manage_session() as session:
            collection_version = session.query(CollectionVersionTable).filter_by(id=collection_version_id.id).one()
            # TODO: alternatively use postgres `array_remove`
            updated_datasets: List[uuid.UUID] = list(collection_version.datasets)
            updated_datasets.remove(uuid.UUID(dataset_version_id.id))
            collection_version.datasets = updated_datasets

    def replace_dataset_in_collection_version(
        self,
        collection_version_id: CollectionVersionId,
        old_dataset_version_id: DatasetVersionId,
        new_dataset_version_id: DatasetVersionId = None,
    ) -> DatasetVersion:
        """
        Replaces an existing mapping between a collection version and a dataset version

        :param collection_version_id: the collection version id
        :param old_dataset_version_id: the dataset version id to be replaced
        :param new_dataset_version_id: the dataset version id to replace with. If None is provide a new
        dataset version will be created.
        """
        # TODO: this method should probably be split into multiple - it contains too much logic
        with self._manage_session() as session:
            collection_id = (
                session.query(CollectionVersionTable.collection_id).filter_by(id=collection_version_id.id).one()[0]
            )  # noqa
            dataset_id = session.query(DatasetVersionTable.dataset_id).filter_by(id=old_dataset_version_id.id).one()[0]
            if new_dataset_version_id is None:
                new_dataset_version_id = DatasetVersionId()
                new_dataset_version = DatasetVersionTable(
                    id=new_dataset_version_id.id,
                    dataset_id=dataset_id,
                    collection_id=collection_id,
                    dataset_metadata=None,
                    artifacts=list(),
                    status=DatasetStatus.empty().to_json(),
                    created_at=datetime.utcnow(),
                )
                session.add(new_dataset_version)
                new_dataset_version = self._hydrate_dataset_version(new_dataset_version)
            else:
                new_dataset_version = self.get_dataset_version(new_dataset_version_id)
                if CollectionId(str(collection_id)) != new_dataset_version.collection_id:
                    raise ValueError(
                        f"Dataset version {new_dataset_version_id} does not belong to collection {collection_id}"
                    )

            collection_version = (
                session.query(CollectionVersionTable).filter_by(id=collection_version_id.id).one()
            )  # noqa
            # This replaces the dataset while preserving the order of datasets
            datasets = list(collection_version.datasets)
            idx = next(i for i, e in enumerate(datasets) if str(e) == old_dataset_version_id.id)
            datasets[idx] = uuid.UUID(new_dataset_version_id.id)
            collection_version.datasets = datasets

            return new_dataset_version

    def get_dataset_mapped_version(
        self, dataset_id: DatasetId, get_tombstoned: bool = False
    ) -> Optional[DatasetVersion]:
        """
        Returns the dataset version mapped to a canonical dataset_id, or None if not existing
        """
        with self._manage_session() as session:
            if get_tombstoned:
                canonical_dataset = session.query(DatasetTable).filter_by(id=dataset_id.id)
            else:
                canonical_dataset = session.query(DatasetTable).filter_by(id=dataset_id.id).filter_by(tombstone=False)
            canonical_dataset = canonical_dataset.one_or_none()

            if canonical_dataset is None:
                return None
            if canonical_dataset.version_id is None:
                return None
            dataset_version = session.query(DatasetVersionTable).filter_by(id=canonical_dataset.version_id).one()
            dataset_version.canonical_dataset = canonical_dataset
            return self._hydrate_dataset_version(dataset_version)

    def get_collection_versions_by_schema(self, schema_version: str, has_wildcards: bool) -> List[CollectionVersion]:
        """
        Returns a list with all collection versions that match the given schema_version. schema_version may contain
         wildcards.
        """
        with self._manage_session() as session:
            if has_wildcards:
                collection_versions = session.query(CollectionVersionTable).filter(
                    CollectionVersionTable.schema_version.like(schema_version)
                )
            else:
                collection_versions = session.query(CollectionVersionTable).filter_by(schema_version=schema_version)
            return [self._row_to_collection_version(row, None) for row in collection_versions.all()]

    def get_previous_dataset_version_id(self, dataset_id: DatasetId) -> Optional[DatasetVersionId]:
        """
        Returns the previously created dataset version for a dataset.
        """
        with self._manage_session() as session:
            version_id = (
                session.query(DatasetVersionTable.id)
                .filter_by(dataset_id=dataset_id.id)
                .order_by(DatasetVersionTable.created_at.desc())
                .offset(1)
                .first()
            )
            if version_id is None:
                return None
            return DatasetVersionId(str(version_id.id))
