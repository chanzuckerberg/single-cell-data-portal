from typing import Iterable, List, Optional, Tuple

from backend.layers.business.entities import (
    CollectionMetadataUpdate,
    CollectionQueryFilter,
    DatasetArtifactDownloadData,
)
from backend.layers.common.entities import (
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
    def get_collections(self, filter: CollectionQueryFilter) -> Iterable[CollectionVersion]:  # type: ignore
        pass

    def get_published_collection_version(self, collection_id: CollectionId) -> Optional[CollectionVersionWithDatasets]:
        pass

    def get_unpublished_collection_version_from_canonical(
        self, collection_id: CollectionId
    ) -> Optional[CollectionVersionWithDatasets]:
        pass

    def get_collection_version(  # type: ignore
        self, version_id: CollectionVersionId, get_tombstoned: bool
    ) -> CollectionVersionWithDatasets:
        pass

    def get_collection_versions_from_canonical(self, collection_id: CollectionId) -> Iterable[CollectionVersion]:  # type: ignore
        pass

    def get_collection_version_from_canonical(
        self, collection_id: CollectionId
    ) -> Optional[CollectionVersionWithDatasets]:
        pass

    def create_collection(  # type: ignore
        self, owner: str, curator_name: str, collection_metadata: CollectionMetadata
    ) -> CollectionVersion:
        pass

    def delete_dataset_versions_from_bucket(self, dataset_version_ids: List[str], bucket: str) -> List[str]:  # type: ignore
        pass

    def delete_all_dataset_versions_from_bucket_for_collection(  # type: ignore
        self, collection_id: CollectionId, bucket: str
    ) -> List[str]:
        pass

    def delete_collection(self, collection_id: CollectionId) -> None:
        pass

    def update_collection_version(self, version_id: CollectionVersionId, body: CollectionMetadataUpdate) -> None:
        pass

    def create_collection_version(self, collection_id: CollectionId) -> CollectionVersion:  # type: ignore
        pass

    def delete_collection_version(self, version_id: CollectionVersionId) -> None:
        pass

    def publish_collection_version(self, version_id: CollectionVersionId) -> None:
        pass

    def ingest_dataset(  # type: ignore
        self,
        collection_version_id: CollectionVersionId,
        url: str,
        file_size: Optional[int],
        existing_dataset_version_id: Optional[DatasetVersionId],
    ) -> Tuple[DatasetVersionId, DatasetId]:
        pass

    def get_all_mapped_datasets(self) -> List[DatasetVersion]:  # type: ignore
        pass

    def get_all_mapped_collection_versions_with_datasets(self) -> List[CollectionVersionWithPublishedDatasets]:  # type: ignore
        pass

    def remove_dataset_version(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId, delete_published: bool
    ) -> None:
        pass

    def set_dataset_metadata(self, dataset_version_id: DatasetVersionId, metadata: DatasetMetadata) -> None:
        pass

    def get_dataset_artifacts(self, dataset_version_id: DatasetVersionId) -> Iterable[DatasetArtifact]:  # type: ignore
        pass

    def get_dataset_artifact_download_data(  # type: ignore
        self, dataset_version_id: DatasetVersionId, artifact_id: DatasetArtifactId
    ) -> DatasetArtifactDownloadData:
        pass

    def update_dataset_version_status(
        self,
        dataset_version_id: DatasetVersionId,
        status_key: DatasetStatusKey,
        new_dataset_status: DatasetStatusGeneric,
        validation_message: Optional[str] = None,
    ) -> None:
        pass

    def add_dataset_artifact(  # type: ignore
        self, dataset_version_id: DatasetVersionId, artifact_type: str, artifact_uri: str
    ) -> DatasetArtifactId:
        pass

    def get_dataset_status(self, dataset_version_id: DatasetVersionId) -> DatasetStatus:  # type: ignore
        pass

    def get_dataset_version(self, dataset_version_id: DatasetVersionId) -> DatasetVersion:  # type: ignore
        pass

    def get_prior_published_versions_for_dataset(self, dataset_id: DatasetId) -> List[PublishedDatasetVersion]:  # type: ignore
        pass

    def get_prior_published_dataset_version(self, dataset_version_id: DatasetVersionId) -> PublishedDatasetVersion:  # type: ignore
        pass

    def get_dataset_version_from_canonical(self, dataset_id: DatasetId, get_tombstoned: bool) -> DatasetVersion:  # type: ignore
        pass

    def get_latest_published_collection_versions_by_schema(  # type: ignore
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
