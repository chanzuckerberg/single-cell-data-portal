from typing import Iterable, List, Optional, Tuple

from backend.layers.business.entities import (
    CollectionMetadataUpdate,
    CollectionQueryFilter,
    DatasetArtifactDownloadData,
    DeprecatedDatasetArtifactDownloadData,
)
from backend.layers.common.entities import (
    CanonicalCollection,
    CollectionId,
    CollectionMetadata,
    CollectionVersion,
    CollectionVersionId,
    CollectionVersionWithDatasets,
    CollectionVersionWithPublishedDatasets,
    DatasetArtifact,
    DatasetArtifactId,
    DatasetId,
    DatasetMetadata,
    DatasetStatus,
    DatasetStatusGeneric,
    DatasetStatusKey,
    DatasetVersion,
    DatasetVersionId,
    PublishedDatasetVersion,
)


class BusinessLogicInterface:
    def get_collections(self, filter: CollectionQueryFilter) -> Iterable[CollectionVersion]:
        pass

    def get_published_collection_version(self, collection_id: CollectionId) -> Optional[CollectionVersionWithDatasets]:
        pass

    def get_unpublished_collection_version_from_canonical(
        self, collection_id: CollectionId
    ) -> Optional[CollectionVersionWithDatasets]:
        pass

    def get_collection_version(
        self, version_id: CollectionVersionId, get_tombstoned: bool
    ) -> CollectionVersionWithDatasets:
        pass

    def get_collection_versions_from_canonical(self, collection_id: CollectionId) -> Iterable[CollectionVersion]:
        pass

    def get_canonical_collection(self, collection_id: CollectionId) -> CanonicalCollection:
        pass

    def get_all_published_collection_versions_from_canonical(
        self, collection_id: CollectionId, get_tombstoned: bool
    ) -> Iterable[CollectionVersionWithDatasets]:
        pass

    def get_collection_version_from_canonical(
        self, collection_id: CollectionId
    ) -> Optional[CollectionVersionWithDatasets]:
        pass

    def create_collection(
        self, owner: str, curator_name: str, collection_metadata: CollectionMetadata
    ) -> CollectionVersion:
        pass

    def delete_artifacts(self, artifacts: List[DatasetArtifact]) -> None:
        pass

    def delete_dataset_versions_from_public_bucket(self, dataset_version_ids: List[str]) -> List[str]:
        pass

    def delete_all_dataset_versions_from_public_bucket_for_collection(self, collection_id: CollectionId) -> List[str]:
        pass

    def get_unpublished_dataset_versions(self, dataset_id: DatasetId) -> List[DatasetVersion]:
        """
        Get all unpublished versions for a Dataset, as determined by the DatasetVersion argument, currently in the
        database. Generally this will be one Dataset version unless a Dataset is in the midst of being replaced, in
        which case there will be two DatasetVersion objects associated with the Dataset. It is also possible that
        through pipeline errors, multiple unpublished DatasetVersionTable rows end up being persisted in the database,
        in which case this function can be used for cleanup.
        """
        pass

    def delete_collection(self, collection_id: CollectionId) -> None:
        pass

    def delete_dataset_version_assets(self, dataset_versions: List[DatasetVersion]) -> None:
        pass

    def update_collection_version(self, version_id: CollectionVersionId, body: CollectionMetadataUpdate) -> None:
        pass

    def create_collection_version(self, collection_id: CollectionId) -> CollectionVersion:
        pass

    def delete_collection_version(self, collection_version: CollectionVersionWithDatasets) -> None:
        pass

    def publish_collection_version(self, version_id: CollectionVersionId) -> None:
        pass

    def ingest_dataset(
        self,
        collection_version_id: CollectionVersionId,
        url: str,
        file_size: Optional[int],
        current_dataset_version_id: Optional[DatasetVersionId],
        start_step_function: bool = False,
    ) -> Tuple[DatasetVersionId, DatasetId]:
        pass

    def get_all_mapped_datasets(self) -> List[DatasetVersion]:
        pass

    def get_all_mapped_collection_versions_with_datasets(self) -> List[CollectionVersionWithPublishedDatasets]:
        pass

    def remove_dataset_version(
        self,
        collection_version_id: CollectionVersionId,
        dataset_version_id: DatasetVersionId,
        delete_published: bool = False,
    ) -> None:
        pass

    def set_dataset_metadata(self, dataset_version_id: DatasetVersionId, metadata: DatasetMetadata) -> None:
        pass

    def get_dataset_artifacts(self, dataset_version_id: DatasetVersionId) -> Iterable[DatasetArtifact]:
        pass

    def get_dataset_artifact_download_data(
        self, dataset_version_id: DatasetVersionId, artifact_id: DatasetArtifactId
    ) -> DatasetArtifactDownloadData:
        pass

    # TODO: Superseded by get_dataset_artifact_download_data. Remove with #5697.
    def get_dataset_artifact_download_data_deprecated(
        self, dataset_version_id: DatasetVersionId, artifact_id: DatasetArtifactId
    ) -> DeprecatedDatasetArtifactDownloadData:
        pass

    def update_dataset_version_status(
        self,
        dataset_version_id: DatasetVersionId,
        status_key: DatasetStatusKey,
        new_dataset_status: DatasetStatusGeneric,
        validation_message: Optional[str] = None,
    ) -> None:
        pass

    def add_dataset_artifact(
        self, dataset_version_id: DatasetVersionId, artifact_type: str, artifact_uri: str
    ) -> DatasetArtifactId:
        pass

    def get_dataset_status(self, dataset_version_id: DatasetVersionId) -> DatasetStatus:
        pass

    def get_dataset_version(self, dataset_version_id: DatasetVersionId) -> DatasetVersion:
        pass

    def get_prior_published_versions_for_dataset(self, dataset_id: DatasetId) -> List[PublishedDatasetVersion]:
        pass

    def get_prior_published_dataset_version(self, dataset_version_id: DatasetVersionId) -> PublishedDatasetVersion:
        pass

    def get_dataset_version_from_canonical(self, dataset_id: DatasetId, get_tombstoned: bool) -> DatasetVersion:
        pass

    def get_latest_published_collection_versions_by_schema(
        self, schema_version: str
    ) -> List[CollectionVersionWithPublishedDatasets]:
        pass

    def restore_previous_dataset_version(
        self, collection_version_id: CollectionVersionId, dataset_id: DatasetId
    ) -> None:
        """
        Restore the previous dataset version for a dataset.
        :param collection_version_id: The collection version to restore the dataset version. It must be in a mutable
        state
        :param dataset_id: The dataset id to restore the previous version of.
        """
