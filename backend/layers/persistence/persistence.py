import json
import logging
import uuid
from contextlib import contextmanager
from datetime import datetime
from typing import Any, Iterable, List, Optional

from sqlalchemy import create_engine
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker

from backend.common.corpora_config import CorporaDbConfig
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
from backend.layers.business.exceptions import CollectionIsPublishedException
from backend.layers.persistence.orm import (
    Collection as CollectionTable,
    CollectionVersion as CollectionVersionTable,
    Dataset as DatasetTable,
    DatasetArtifact as DatasetArtifactTable,
    DatasetVersion as DatasetVersionTable,
)
from backend.layers.persistence.persistence_interface import DatabaseProviderInterface, PersistenceException

logger = logging.getLogger(__name__)


class DatabaseProvider(DatabaseProviderInterface):
    def __init__(self, database_uri: str = None, schema_name: str = "persistence_schema") -> None:
        if not database_uri:
            database_uri = CorporaDbConfig().database_uri
        self._engine = create_engine(database_uri, connect_args={"connect_timeout": 5})
        self._session_maker = sessionmaker(bind=self._engine)
        try:
            self._create_schema(schema_name)
        except Exception:
            pass

    def _drop_schema(self, schema_name: str):
        from sqlalchemy.schema import DropSchema

        self._engine.execute(DropSchema(schema_name, cascade=True))

    def _create_schema(self, schema_name: str):
        from sqlalchemy.schema import CreateSchema
        from backend.layers.persistence.orm import metadata

        self._engine.execute(CreateSchema(schema_name))
        metadata.schema = schema_name
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
            raise PersistenceException("Failed to commit.")
        finally:
            if session is not None:
                session.close()

    @staticmethod
    def _generate_id():
        return str(uuid.uuid4())

    def _row_to_collection_version(self, row: Any, canonical_collection: CanonicalCollection) -> CollectionVersion:
        return CollectionVersion(
            collection_id=CollectionId(str(row.collection_id)),
            version_id=CollectionVersionId(str(row.version_id)),
            owner=row.owner,
            curator_name=row.curator_name,
            metadata=CollectionMetadata.from_json(row.metadata),
            publisher_metadata=None if row.publisher_metadata is None else json.loads(row.publisher_metadata),
            datasets=[DatasetVersionId(str(id)) for id in row.datasets],
            published_at=row.published_at,
            created_at=row.created_at,
            canonical_collection=canonical_collection,
        )

    def _row_to_collection_version_with_datasets(
        self, row: Any, canonical_collection: CanonicalCollection, datasets: List[DatasetVersion]
    ) -> CollectionVersionWithDatasets:
        return CollectionVersionWithDatasets(
            collection_id=CollectionId(str(row.collection_id)),
            version_id=CollectionVersionId(str(row.version_id)),
            owner=row.owner,
            curator_name=row.curator_name,
            metadata=CollectionMetadata.from_json(row.metadata),
            publisher_metadata=None if row.publisher_metadata is None else json.loads(row.publisher_metadata),
            datasets=datasets,
            published_at=row.published_at,
            created_at=row.created_at,
            canonical_collection=canonical_collection,
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
            artifacts,
            row.created_at,
            canonical_dataset,
        )

    def _hydrate_dataset_version(self, dataset_version: DatasetVersionTable) -> DatasetVersion:
        """
        Populates canonical_dataset, artifacts, and status for DatasetVersionRow
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
                collection.tombstoned,
            )

    def get_canonical_dataset(self, dataset_id: DatasetId) -> CanonicalDataset:
        with self._manage_session() as session:
            dataset = session.query(DatasetTable).filter_by(dataset_id=dataset_id.id).one_or_none()
            if dataset is None:
                return None
            return CanonicalDataset(dataset_id, DatasetVersionId(str(dataset.dataset_version_id)), dataset.published_at)

    def create_canonical_collection(
        self, owner: str, curator_name: str, collection_metadata: CollectionMetadata
    ) -> CollectionVersion:
        """
        Creates a new canonical collection, generating a canonical collection_id and a new version_id.
        Returns the newly created CollectionVersion
        """
        collection_id = CollectionId((self._generate_id()))
        version_id = CollectionVersionId((self._generate_id()))
        now = datetime.utcnow()
        canonical_collection = CollectionTable(
            id=collection_id.id,
            version_id=None,
            tombstoned=False,
            originally_published_at=None,
        )

        collection_version_row = CollectionVersionTable(
            collection_id=collection_id.id,
            version_id=version_id.id,
            owner=owner,
            curator_name=curator_name,
            metadata=collection_metadata.to_json(),
            publisher_metadata=None,
            published_at=None,
            created_at=now,
            datasets=list(),
        )

        with self._manage_session() as session:
            session.add(canonical_collection)
            session.add(collection_version_row)

            return self._row_to_collection_version(
                collection_version_row, CanonicalCollection(collection_id, None, None, False)
            )

    def get_collection_version(self, version_id: CollectionVersionId) -> CollectionVersion:
        """
        Retrieves a specific collection version by id
        """
        with self._manage_session() as session:
            collection_version = session.query(CollectionVersionTable).filter_by(version_id=version_id.id).one_or_none()
            if collection_version is None:
                return None
            collection_id = CollectionId(str(collection_version.collection_id))
            canonical_collection = self.get_canonical_collection(collection_id)
            return self._row_to_collection_version(collection_version, canonical_collection)

    def _get_datasets(self, ids: List[DatasetVersionId]):
        # TODO: can be optimized via in queries, if necessary
        return [self.get_dataset_version(id) for id in ids]

    def get_collection_version_with_datasets(self, version_id: CollectionVersionId) -> CollectionVersionWithDatasets:
        """
        Retrieves a specific collection version by id, with datasets
        """
        with self._manage_session() as session:
            collection_version = session.query(CollectionVersionTable).filter_by(version_id=version_id.id).one_or_none()
            if collection_version is None:
                return None
            collection_id = CollectionId(str(collection_version.collection_id))
            canonical_collection = self.get_canonical_collection(collection_id)
            datasets = self._get_datasets([DatasetVersionId(str(id)) for id in collection_version.datasets])
            return self._row_to_collection_version_with_datasets(collection_version, canonical_collection, datasets)

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
            collection_version = session.query(CollectionVersionTable).filter_by(version_id=version_id).one()
            canonical_collection = self.get_canonical_collection(collection_id)
            datasets = self._get_datasets([DatasetVersionId(str(id)) for id in collection_version.datasets])
            return self._row_to_collection_version_with_datasets(collection_version, canonical_collection, datasets)

    def get_all_versions_for_collection(self, collection_id: CollectionId) -> List[CollectionVersionWithDatasets]:
        """
        Retrieves all versions for a specific collections, without filtering
        """
        with self._manage_session() as session:
            version_rows = session.query(CollectionVersionTable).filter_by(collection_id=collection_id.id).all()
            canonical_collection = self.get_canonical_collection(collection_id)
            versions = list()
            for i in range(len(version_rows)):
                datasets = self._get_datasets([DatasetVersionId(str(id)) for id in version_rows[i].datasets])
                version = self._row_to_collection_version_with_datasets(version_rows[i], canonical_collection, datasets)
                versions.append(version)
            return versions

    def get_all_collections_versions(self) -> Iterable[CollectionVersion]:
        """
        Retrieves all versions of all collections.
        TODO: for performance reasons, it might be necessary to add a filtering parameter here.
        """
        with self._manage_session() as session:
            versions = session.query(CollectionVersionTable).all()
            # TODO: do we need to hydrate versions with canonical collections? would require a join or many lookup calls
            return [self._row_to_collection_version(v, None) for v in versions]

    def get_all_mapped_collection_versions(self, get_tombstoned: bool = False) -> Iterable[CollectionVersion]:
        """
        Retrieves all the collection versions that are mapped to a canonical collection.
        """
        with self._manage_session() as session:
            if get_tombstoned:
                canonical_collections = (
                    session.query(CollectionTable).filter(CollectionTable.version_id.isnot(None)).all()
                )  # noqa
            else:
                canonical_collections = (
                    session.query(CollectionTable)
                    .filter(CollectionTable.version_id.isnot(None))
                    .filter_by(tombstoned=False)
                    .all()
                )

            mapped_version_ids = [i.version_id for i in canonical_collections]
            versions = (
                session.query(CollectionVersionTable)
                .filter(CollectionVersionTable.version_id.in_(mapped_version_ids))
                .all()
            )  # noqa

            for version in versions:
                # TODO: should be optimized using a map
                canonical_row = next(cc for cc in canonical_collections if cc.version_id == version.version_id)
                canonical = CanonicalCollection(
                    CollectionId(str(canonical_row.id)),
                    CollectionVersionId(str(canonical_row.version_id)),
                    canonical_row.originally_published_at,
                    canonical_row.tombstoned,
                )
                print("CCCCCCCCCCCCC", canonical)

                yield self._row_to_collection_version(version, canonical)

    def delete_canonical_collection(self, collection_id: CollectionId) -> None:
        """
        Deletes (tombstones) a canonical collection.
        """
        with self._manage_session() as session:
            canonical_collection = session.query(CollectionTable).filter_by(id=collection_id.id).one_or_none()
            if canonical_collection:
                canonical_collection.tombstoned = True

    def save_collection_metadata(
        self, version_id: CollectionVersionId, collection_metadata: CollectionMetadata
    ) -> None:
        """
        Saves collection metadata for a collection version
        """
        with self._manage_session() as session:
            version = session.query(CollectionVersionTable).filter_by(version_id=version_id.id).one()
            version.metadata = collection_metadata.to_json()

    def save_collection_publisher_metadata(
        self, version_id: CollectionVersionId, publisher_metadata: Optional[dict]
    ) -> None:
        """
        Saves publisher metadata for a collection version. Specify None to remove it
        """
        with self._manage_session() as session:
            version = session.query(CollectionVersionTable).filter_by(version_id=version_id.id).one()
            version.publisher_metadata = json.dumps(publisher_metadata)

    def add_collection_version(self, collection_id: CollectionId) -> CollectionVersionId:
        """
        Adds a collection version to an existing canonical collection. The new version copies all data from
         the previous version except version_id and datetime-based fields (i.e. created_at, published_at)
        Returns the new version id.
        """
        with self._manage_session() as session:
            current_version_id = session.query(CollectionTable.version_id).filter_by(id=collection_id.id).one()[0]
            current_version = session.query(CollectionVersionTable).filter_by(version_id=current_version_id).one()
            new_version_id = self._generate_id()
            new_version = CollectionVersionTable(
                version_id=new_version_id,
                collection_id=collection_id.id,
                metadata=current_version.metadata,
                owner=current_version.owner,
                curator_name=current_version.curator_name,
                publisher_metadata=current_version.publisher_metadata,
                published_at=None,
                created_at=datetime.utcnow(),
                datasets=current_version.datasets,
            )
            session.add(new_version)
            return CollectionVersionId(new_version_id)

    def delete_collection_version(self, version_id: CollectionVersionId) -> None:
        """
        Deletes a collection version, if it is unpublished.
        """
        with self._manage_session() as session:
            version = session.query(CollectionVersionTable).filter_by(version_id=version_id.id).one_or_none()
            if version:
                if version.published_at:
                    raise CollectionIsPublishedException(f"Published Collection Version {version_id} cannot be deleted")
                session.delete(version)

    def finalize_collection_version(
        self,
        collection_id: CollectionId,
        version_id: CollectionVersionId,
        published_at: Optional[datetime] = None,
    ) -> None:
        """
        Finalizes a collection version. This is equivalent to calling:
        1. update_collection_version_mapping
        2. set_collection_version_published_at
        3. finalize_dataset_versions
        """
        published_at = published_at if published_at else datetime.utcnow()
        self.update_collection_version_mapping(collection_id, version_id, published_at)
        self.set_collection_version_published_at(version_id, published_at)
        self.finalize_dataset_versions(version_id, published_at)

    def update_collection_version_mapping(
        self, collection_id: CollectionId, version_id: CollectionVersionId, published_at: datetime
    ) -> None:
        """
        Updates the mapping between the canonical collection `collection_id` and its `version_id`
        """
        with self._manage_session() as session:
            collection = session.query(CollectionTable).filter_by(id=collection_id.id).one()
            collection.version_id = version_id.id
            if collection.originally_published_at is None:
                collection.originally_published_at = published_at

    def set_collection_version_published_at(self, version_id: CollectionVersionId, published_at: datetime) -> None:
        """
        Sets the `published_at` datetime for a collection version
        """
        with self._manage_session() as session:
            collection_version = session.query(CollectionVersionTable).filter_by(version_id=version_id.id).one()
            collection_version.published_at = published_at

    def finalize_dataset_versions(self, version_id: CollectionVersionId, published_at: datetime) -> None:
        """
        1. Updates the mapping between the canonical dataset and the latest published dataset version for each dataset
        in a CollectionVersion
        2. Sets the each canonical dataset's 'published_at' if not previously set
        """
        with self._manage_session() as session:
            dataset_version_ids = (
                session.query(CollectionVersionTable.datasets).filter_by(version_id=version_id.id).one()[0]
            )
            for dataset_version, dataset in (
                session.query(DatasetVersionTable, DatasetTable)
                .filter(DatasetVersionTable.dataset_id == DatasetTable.dataset_id)
                .filter(DatasetVersionTable.version_id.in_(dataset_version_ids))
                .all()
            ):
                dataset.version_id = dataset_version.version_id
                if dataset.published_at is None:
                    dataset.published_at = published_at

    def get_dataset_version(self, dataset_version_id: DatasetVersionId) -> DatasetVersion:
        """
        Returns a dataset version by id.
        """
        with self._manage_session() as session:
            dataset_version = (
                session.query(DatasetVersionTable).filter_by(version_id=dataset_version_id.id).one_or_none()
            )
            if dataset_version is None:
                return None
            return self._hydrate_dataset_version(dataset_version)

    def get_all_versions_for_dataset(self, dataset_id: DatasetId) -> List[DatasetVersion]:
        """
        Returns all dataset versions for a canonical dataset_id
        """
        dataset = self.get_canonical_dataset(dataset_id)
        with self._manage_session() as session:
            dataset_versions = session.query(DatasetVersionTable).filter_by(dataset_id=dataset_id.id).all()
            for i in range(len(dataset_versions)):
                dataset_versions[i] = self._row_to_dataset_version(dataset_versions[i], dataset)
            return dataset_versions

    def get_all_datasets(self) -> Iterable[DatasetVersion]:
        """
        Returns all dataset versions.
        # TODO: Add filtering (tombstoned? remove orphaned datasets? canonical only? published?)
        """
        active_collections = self.get_all_mapped_collection_versions()
        active_datasets = [i.id for s in [c.datasets for c in active_collections] for i in s]

        # TODO: this is doing N fetches - optimize this using INs or by loading the other tables in memory
        acc = []
        with self._manage_session() as session:
            for version in session.query(DatasetVersionTable).all():  # noqa
                if str(version.version_id) in active_datasets:
                    acc.append(self._hydrate_dataset_version(version))
        return acc

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
            print("zzz", artifact_ids)
        return self.get_dataset_artifacts(artifact_ids[0])

    def create_canonical_dataset(self, collection_version_id: CollectionVersionId) -> DatasetVersion:
        """
        Initializes a canonical dataset, generating a dataset_id and a dataset_version_id.
        Returns the newly created DatasetVersion.
        """
        with self._manage_session() as session:
            collection_id = (
                session.query(CollectionVersionTable.collection_id)
                .filter_by(version_id=collection_version_id.id)
                .one()[0]
            )
        dataset_id = DatasetId(self._generate_id())
        dataset_version_id = DatasetVersionId(self._generate_id())
        canonical_dataset = DatasetTable(
            dataset_id=dataset_id.id, dataset_version_id=dataset_version_id.id, published_at=None
        )
        dataset_version = DatasetVersionTable(
            version_id=dataset_version_id.id,
            dataset_id=dataset_id.id,
            collection_id=collection_id,
            metadata=None,
            artifacts=list(),
            status=DatasetStatus.empty().to_json(),
            created_at=datetime.utcnow(),
        )

        with self._manage_session() as session:
            session.add(canonical_dataset)
            session.add(dataset_version)
            return self._row_to_dataset_version(
                dataset_version, CanonicalDataset(dataset_id, dataset_version_id, None), []
            )

    def add_dataset_artifact(
        self, version_id: DatasetVersionId, artifact_type: DatasetArtifactType, artifact_uri: str
    ) -> DatasetArtifactId:
        """
        Adds a dataset artifact to an existing dataset version.
        """
        artifact_id = DatasetArtifactId(self._generate_id())
        artifact = DatasetArtifactTable(id=artifact_id.id, type=artifact_type, uri=artifact_uri)
        with self._manage_session() as session:
            session.add(artifact)
            dataset_version = session.query(DatasetVersionTable).filter_by(version_id=version_id.id).one()
            artifacts = list(dataset_version.artifacts)
            artifacts.append(uuid.UUID(artifact_id.id))
            dataset_version.artifacts = artifacts
        return artifact_id

    def update_dataset_processing_status(self, version_id: DatasetVersionId, status: DatasetProcessingStatus) -> None:
        """
        Updates the processing status for a dataset version.
        """
        with self._manage_session() as session:
            dataset_version = session.query(DatasetVersionTable).filter_by(version_id=version_id.id).one()
            dataset_version_status = json.loads(dataset_version.status)
            dataset_version_status["processing_status"] = status.value
            dataset_version.status = json.dumps(dataset_version_status)

    def update_dataset_validation_status(self, version_id: DatasetVersionId, status: DatasetValidationStatus) -> None:
        """
        Updates the validation status for a dataset version.
        """
        with self._manage_session() as session:
            dataset_version = session.query(DatasetVersionTable).filter_by(version_id=version_id.id).one()
            dataset_version_status = json.loads(dataset_version.status)
            dataset_version_status["validation_status"] = status.value
            dataset_version.status = json.dumps(dataset_version_status)

    def update_dataset_upload_status(self, version_id: DatasetVersionId, status: DatasetUploadStatus) -> None:
        """
        Updates the upload status for a dataset version.
        """
        with self._manage_session() as session:
            dataset_version = session.query(DatasetVersionTable).filter_by(version_id=version_id.id).one()
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
            dataset_version = session.query(DatasetVersionTable).filter_by(version_id=version_id.id).one()
            dataset_version_status = json.loads(dataset_version.status)
            dataset_version_status[status_type] = status.value
            dataset_version.status = json.dumps(dataset_version_status)

    def update_dataset_validation_message(self, version_id: DatasetVersionId, validation_message: str) -> None:
        with self._manage_session() as session:
            dataset_version = session.query(DatasetVersionTable).filter_by(version_id=version_id.id).one()
            dataset_version_status = json.loads(dataset_version.status)
            dataset_version_status["validation_message"] = validation_message
            dataset_version.status = json.dumps(dataset_version_status)

    def get_dataset_version_status(self, version_id: DatasetVersionId) -> DatasetStatus:
        """
        Returns the status for a dataset version
        """
        with self._manage_session() as session:
            status = session.query(DatasetVersionTable.status).filter_by(version_id=version_id.id).one()
        return DatasetStatus.from_json(status[0])

    def set_dataset_metadata(self, version_id: DatasetVersionId, metadata: DatasetMetadata) -> None:
        """
        Sets the metadata for a dataset version
        """
        with self._manage_session() as session:
            dataset_version = session.query(DatasetVersionTable).filter_by(version_id=version_id.id).one()
            dataset_version.metadata = metadata.to_json()

    def add_dataset_to_collection_version_mapping(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId
    ) -> None:
        """
        Adds a mapping between an existing collection version and a dataset version
        """
        with self._manage_session() as session:
            collection_version = (
                session.query(CollectionVersionTable).filter_by(version_id=collection_version_id.id).one()
            )
            # TODO: alternatively use postgres `array_append`
            # TODO: make sure that the UUID conversion works
            updated_datasets = list(collection_version.datasets)
            # print("before", updated_datasets)
            updated_datasets.append(uuid.UUID(dataset_version_id.id))
            collection_version.datasets = updated_datasets
            # print("after", updated_datasets)

    def delete_dataset_from_collection_version(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId
    ) -> None:
        """
        Removes a mapping between a collection version and a dataset version
        """
        with self._manage_session() as session:
            collection_version = (
                session.query(CollectionVersionTable).filter_by(version_id=collection_version_id.id).one()
            )
            # TODO: alternatively use postgres `array_remove`
            updated_datasets = list(collection_version.datasets)
            updated_datasets.remove(uuid.UUID(dataset_version_id.id))
            collection_version.datasets = updated_datasets

    def replace_dataset_in_collection_version(
        self, collection_version_id: CollectionVersionId, old_dataset_version_id: DatasetVersionId
    ) -> DatasetVersion:
        """
        Replaces an existing mapping between a collection version and a dataset version
        """
        # TODO: this method should probably be split into multiple - it contains too much logic
        with self._manage_session() as session:
            collection_id = (
                session.query(CollectionVersionTable.collection_id)
                .filter_by(version_id=collection_version_id.id)
                .one()[0]
            )  # noqa
            dataset_id = (
                session.query(DatasetVersionTable.dataset_id).filter_by(version_id=old_dataset_version_id.id).one()[0]
            )
            new_dataset_version_id = DatasetVersionId(self._generate_id())
            new_dataset_version = DatasetVersionTable(
                version_id=new_dataset_version_id.id,
                dataset_id=dataset_id,
                collection_id=collection_id,
                metadata=None,
                artifacts=list(),
                status=DatasetStatus.empty().to_json(),
                created_at=datetime.utcnow(),
            )
            session.add(new_dataset_version)

            collection_version = (
                session.query(CollectionVersionTable).filter_by(version_id=collection_version_id.id).one()
            )  # noqa
            # This replaces the dataset while preserving the order of datasets
            datasets = list(collection_version.datasets)
            idx = next(i for i, e in enumerate(datasets) if str(e) == old_dataset_version_id.id)
            datasets[idx] = uuid.UUID(new_dataset_version_id.id)
            collection_version.datasets = datasets

            return self._hydrate_dataset_version(new_dataset_version)

    def get_dataset_mapped_version(self, dataset_id: DatasetId) -> Optional[DatasetVersion]:
        """
        Returns the dataset version mapped to a canonical dataset_id, or None if not existing
        """
        with self._manage_session() as session:
            canonical_dataset = session.query(DatasetTable).filter_by(dataset_id=dataset_id.id).one_or_none()
            if canonical_dataset is None:
                return None
            if canonical_dataset.dataset_version_id is None:
                return None
            dataset_version = (
                session.query(DatasetVersionTable).filter_by(version_id=canonical_dataset.dataset_version_id).one()
            )
            dataset_version.canonical_dataset = canonical_dataset
            return self._hydrate_dataset_version(dataset_version)
