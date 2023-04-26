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

    def get_collection_version(self, version_id: CollectionVersionId) -> CollectionVersionWithDatasets:
        pass

    def get_collection_versions_from_canonical(self, collection_id: CollectionId) -> Iterable[CollectionVersion]:
        pass

    def get_collection_version_from_canonical(
        self, collection_id: CollectionId
    ) -> Optional[CollectionVersionWithDatasets]:
        pass

    def create_collection(
        self, owner: str, curator_name: str, collection_metadata: CollectionMetadata
    ) -> CollectionVersion:
        pass

    def delete_datasets_from_bucket(self, collection_id: CollectionId, bucket: str) -> List[str]:
        pass

    def delete_collection(self, collection_id: CollectionId) -> None:
        pass

    def update_collection_version(self, version_id: CollectionVersionId, body: CollectionMetadataUpdate) -> None:
        pass

    def create_collection_version(self, collection_id: CollectionId) -> CollectionVersion:
        pass

    def delete_collection_version(self, version_id: CollectionVersionId) -> None:
        pass

    def publish_collection_version(self, version_id: CollectionVersionId) -> None:
        pass

    def ingest_dataset(
        self,
        collection_version_id: CollectionVersionId,
        url: str,
        file_size: Optional[int],
        existing_dataset_version_id: Optional[DatasetVersionId],
    ) -> Tuple[DatasetVersionId, DatasetId]:
        pass

    def get_all_mapped_datasets(self) -> List[DatasetVersion]:
        pass

    def get_all_mapped_collection_versions_with_datasets(self) -> List[CollectionVersionWithPublishedDatasets]:
        pass

    def remove_dataset_version(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId
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

    def get_dataset_version_from_canonical(self, dataset_id: DatasetId) -> DatasetVersion:
        pass
