import copy
import logging
import os
from collections import defaultdict
from datetime import datetime
from functools import reduce
from typing import Dict, Iterable, List, Optional, Set, Tuple

from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.business.entities import (
    CollectionMetadataUpdate,
    CollectionQueryFilter,
    DatasetArtifactDownloadData,
)
from backend.layers.business.exceptions import (
    ArtifactNotFoundException,
    CollectionCreationException,
    CollectionDeleteException,
    CollectionIsPublicException,
    CollectionIsPublishedException,
    CollectionNotFoundException,
    CollectionPublishException,
    CollectionUpdateException,
    CollectionVersionException,
    DatasetIngestException,
    DatasetInWrongStatusException,
    DatasetIsPrivateException,
    DatasetIsTombstonedException,
    DatasetNotFoundException,
    DatasetUpdateException,
    DatasetVersionNotFoundException,
    InvalidURIException,
    MaxFileSizeExceededException,
    NoPreviousDatasetVersionException,
)
from backend.layers.common import validation
from backend.layers.common.cleanup import sanitize
from backend.layers.common.entities import (
    CanonicalCollection,
    CollectionId,
    CollectionLinkType,
    CollectionMetadata,
    CollectionVersion,
    CollectionVersionId,
    CollectionVersionWithDatasets,
    CollectionVersionWithPublishedDatasets,
    DatasetArtifact,
    DatasetArtifactId,
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetId,
    DatasetMetadata,
    DatasetProcessingStatus,
    DatasetStatus,
    DatasetStatusGeneric,
    DatasetStatusKey,
    DatasetUploadStatus,
    DatasetValidationStatus,
    DatasetVersion,
    DatasetVersionId,
    Link,
    PublishedDatasetVersion,
)
from backend.layers.common.helpers import (
    get_published_at_and_collection_version_id_else_not_found,
)
from backend.layers.common.regex import S3_URI_REGEX
from backend.layers.persistence.persistence_interface import DatabaseProviderInterface
from backend.layers.thirdparty.crossref_provider import (
    CrossrefDOINotFoundException,
    CrossrefException,
    CrossrefProviderInterface,
)
from backend.layers.thirdparty.s3_exceptions import S3DeleteException
from backend.layers.thirdparty.s3_provider_interface import S3ProviderInterface
from backend.layers.thirdparty.step_function_provider import StepFunctionProviderInterface
from backend.layers.thirdparty.uri_provider import UriProviderInterface

logger = logging.getLogger(__name__)


class BusinessLogic(BusinessLogicInterface):
    database_provider: DatabaseProviderInterface
    crossref_provider: CrossrefProviderInterface
    step_function_provider: StepFunctionProviderInterface
    s3_provider: S3ProviderInterface
    uri_provider: UriProviderInterface

    def __init__(
        self,
        database_provider: DatabaseProviderInterface,
        crossref_provider: CrossrefProviderInterface,
        step_function_provider: StepFunctionProviderInterface,
        s3_provider: S3ProviderInterface,
        uri_provider: UriProviderInterface,
    ) -> None:
        self.crossref_provider = crossref_provider
        self.database_provider = database_provider
        self.step_function_provider = step_function_provider
        self.s3_provider = s3_provider
        self.uri_provider = uri_provider
        super().__init__()

    def _get_publisher_metadata(self, doi: str, errors: list) -> Optional[dict]:
        """
        Retrieves publisher metadata from Crossref.
        """
        try:
            return self.crossref_provider.fetch_metadata(doi)
        except CrossrefDOINotFoundException:
            errors.append({"link_type": CollectionLinkType.DOI, "reason": "DOI cannot be found on Crossref"})
        except CrossrefException as e:
            logging.warning(f"CrossrefException on create_collection: {e}. Will ignore metadata.")
            return None

    def create_collection(
        self, owner: str, curator_name: str, collection_metadata: CollectionMetadata
    ) -> CollectionVersion:
        """
        Creates a collection using the specified metadata. If a DOI is defined, will also
        retrieve publisher metadata from Crossref and add it to the collection.
        """

        sanitize(collection_metadata)

        # Check metadata is valid
        errors = []
        validation.verify_collection_metadata(collection_metadata, errors)

        # TODO: Maybe switch link.type to be an enum
        doi = next((link.uri for link in collection_metadata.links if link.type == "DOI"), None)

        publisher_metadata = self._get_publisher_metadata(doi, errors) if doi is not None else None

        if errors:
            raise CollectionCreationException(errors)

        created_version = self.database_provider.create_canonical_collection(owner, curator_name, collection_metadata)

        # TODO: can collapse with `create_canonical_collection`
        if publisher_metadata:
            self.database_provider.save_collection_publisher_metadata(created_version.version_id, publisher_metadata)
            # Add this to the returned object, to save another `get_collection` call
            created_version.publisher_metadata = publisher_metadata

        return created_version

    def get_published_collection_version(self, collection_id: CollectionId) -> Optional[CollectionVersionWithDatasets]:
        """
        Returns the published collection version that belongs to a canonical collection.
        Returns None if no published collection exists
        """
        return self.database_provider.get_collection_mapped_version(collection_id)

    def get_latest_published_collection_versions_by_schema(
        self, schema_version: str
    ) -> List[CollectionVersionWithPublishedDatasets]:
        """
        Returns a list with the latest published collection version that matches the given schema_version, for each
        canonical collection
        """
        has_wildcards = "_" in schema_version
        collection_versions = self.database_provider.get_collection_versions_by_schema(schema_version, has_wildcards)

        # for each published canonical collection, map its most recently published collection version
        collections = dict()
        for collection_version in collection_versions:
            if collection_version.published_at is None:
                continue
            canonical_collection_id = collection_version.collection_id.id
            if canonical_collection_id not in collections:
                collections[canonical_collection_id] = collection_version
            else:
                if collection_version.published_at > collections[canonical_collection_id].published_at:
                    collections[canonical_collection_id] = collection_version

        # for each mapped collection version, populate its dataset versions' details
        latest_collection_versions = list(collections.values())
        dataset_version_ids = [
            dataset_version_id
            for collection in latest_collection_versions
            for dataset_version_id in collection.datasets
        ]
        dataset_versions = self.database_provider.get_dataset_versions_by_id(dataset_version_ids)
        return self._map_collection_version_to_published_dataset_versions(latest_collection_versions, dataset_versions)

    def get_unpublished_collection_version_from_canonical(
        self, collection_id: CollectionId
    ) -> Optional[CollectionVersionWithDatasets]:
        """
        Given a canonical collection_id, retrieves its latest unpublished version
        """
        latest = datetime.fromtimestamp(0)
        unpublished_collection = None
        for collection in self.get_collection_versions_from_canonical(collection_id):
            if collection.published_at is None and collection.created_at > latest:
                latest = collection.created_at
                unpublished_collection = collection
        return unpublished_collection

    def get_collection_version(
        self, version_id: CollectionVersionId, get_tombstoned: bool = False
    ) -> CollectionVersionWithDatasets:
        """
        Returns a specific collection version by id
        """
        return self.database_provider.get_collection_version_with_datasets(version_id, get_tombstoned=get_tombstoned)

    def get_collection_versions_from_canonical(
        self, collection_id: CollectionId
    ) -> Iterable[CollectionVersionWithDatasets]:
        """
        Returns all the collection versions connected to a canonical collection
        """
        return self.database_provider.get_all_versions_for_collection(collection_id)

    def get_canonical_collection(self, collection_id: CollectionId) -> CanonicalCollection:
        return self.database_provider.get_canonical_collection(collection_id)

    def get_all_published_collection_versions_from_canonical(
        self, collection_id: CollectionId, get_tombstoned: bool = False
    ) -> Iterable[CollectionVersionWithDatasets]:
        """
        Returns all *published* collection versions for a canonical collection
        """
        all_versions = self.database_provider.get_all_versions_for_collection(
            collection_id, get_tombstoned=get_tombstoned
        )
        for c_version in all_versions:
            if c_version.published_at:
                yield c_version

    def get_collection_version_from_canonical(
        self, collection_id: CollectionId, get_tombstoned: bool = False
    ) -> Optional[CollectionVersionWithDatasets]:
        """
        Returns the published collection version mapped to a canonical collection, if available.
        Otherwise will return the active unpublished version.
        """
        version = self.get_published_collection_version(collection_id)
        if version is None:
            version = self.get_unpublished_collection_version_from_canonical(collection_id)
        return version

    def get_collections(self, filter: CollectionQueryFilter) -> Iterable[CollectionVersion]:
        """
        Returns an iterable with all the collections matching `filter`
        """

        # TODO: instead of `is_published`, we should probably call this `is_active_and_published`
        if filter.is_published is True:
            iterable = self.database_provider.get_all_mapped_collection_versions()
        else:
            iterable = self.database_provider.get_all_collections_versions()

        def predicate(version: CollectionVersion):
            # Combines all the filters into a single predicate, so that filtering can happen in a single iteration
            published = (
                filter.is_published is None
                or (filter.is_published is True and version.published_at is not None)
                or (filter.is_published is False and version.published_at is None)
            )
            owner = filter.owner is None or filter.owner == version.owner
            curator = filter.curator_name is None or filter.curator_name == version.curator_name
            return published and owner and curator

        for collection_version in iterable:
            if predicate(collection_version):
                yield collection_version

    def update_collection_version(
        self, version_id: CollectionVersionId, body: CollectionMetadataUpdate, ignore_doi_update: bool = False
    ) -> None:
        """
        Updates a collection version by replacing parts of its metadata.
        If the DOI in the links changed, it should also update its publisher metadata.
        If `ignore_doi_update` is set to True, no DOI updates should be issued

        """

        # TODO: CollectionMetadataUpdate should probably be used for collection creation as well
        # TODO: link.type should DEFINITELY move to an enum. pylance will help with the refactor
        sanitize(body)
        errors = []

        # Check metadata
        validation.verify_collection_metadata_update(body, errors)

        current_version = self.get_collection_version(version_id)
        if current_version.published_at is not None:
            raise CollectionUpdateException(["Cannot update a published collection"])

        # Determine if publisher metadata should be unset, ignored or set at the end of the method.
        # Note: the update needs to be done at the end to ensure atomicity
        unset_publisher_metadata = False
        publisher_metadata_to_set = None

        if not ignore_doi_update:
            # Determine if the DOI has changed
            old_doi = next((link.uri for link in current_version.metadata.links if link.type == "DOI"), None)
            new_doi = (
                None if body.links is None else next((link.uri for link in body.links if link.type == "DOI"), None)
            )

            if old_doi and new_doi is None:
                # If the DOI was deleted, remove the publisher_metadata field
                unset_publisher_metadata = True
            elif (new_doi is not None) and new_doi != old_doi:
                # If the DOI has changed, fetch and update the metadata
                publisher_metadata_to_set = self._get_publisher_metadata(new_doi, errors)

        if errors:
            raise CollectionUpdateException(errors)

        # Merge the updated fields in the existing metadata object. Use a copy to ensure immutability.
        new_metadata = copy.deepcopy(current_version.metadata)
        for field in vars(body):
            if hasattr(body, field) and (value := getattr(body, field)) is not None:
                if isinstance(value, str):
                    value.strip()
                if isinstance(value, Link):
                    value.strip_fields()
                setattr(new_metadata, field, value)

        # Issue all updates
        if unset_publisher_metadata:
            self.database_provider.save_collection_publisher_metadata(version_id, None)
        elif publisher_metadata_to_set is not None:
            self.database_provider.save_collection_publisher_metadata(version_id, publisher_metadata_to_set)
        self.database_provider.save_collection_metadata(version_id, new_metadata)

    def _assert_collection_version_unpublished(
        self, collection_version_id: CollectionVersionId
    ) -> CollectionVersionWithDatasets:
        """
        Ensures a collection version exists and is unpublished.
        This method should be called every time an update to a collection version is requested,
        since published collection versions are not allowed any changes.
        """
        collection_version = self.database_provider.get_collection_version_with_datasets(collection_version_id)
        if collection_version is None:
            raise CollectionNotFoundException([f"Collection version {collection_version_id.id} does not exist"])
        if collection_version.published_at is not None:
            raise CollectionIsPublishedException([f"Collection version {collection_version_id.id} is published"])
        return collection_version

    def create_empty_dataset(self, collection_version_id: CollectionVersionId) -> Tuple[DatasetVersionId, DatasetId]:
        """
        Creates an empty dataset that can be later used for ingestion
        """
        self._assert_collection_version_unpublished(collection_version_id)

        new_dataset_version = self.database_provider.create_canonical_dataset(collection_version_id)
        # adds new dataset version to collection version
        self.database_provider.add_dataset_to_collection_version_mapping(
            collection_version_id, new_dataset_version.version_id
        )

        self.database_provider.update_dataset_upload_status(new_dataset_version.version_id, DatasetUploadStatus.NA)
        self.database_provider.update_dataset_processing_status(
            new_dataset_version.version_id, DatasetProcessingStatus.INITIALIZED
        )

        return (new_dataset_version.version_id, new_dataset_version.dataset_id)

    # TODO: Alternatives: 1) return DatasetVersion 2) Return a new class
    def ingest_dataset(
        self,
        collection_version_id: CollectionVersionId,
        url: str,
        file_size: Optional[int],
        existing_dataset_version_id: Optional[DatasetVersionId],
        start_step_function: bool = True,
    ) -> Tuple[DatasetVersionId, DatasetId]:
        """
        Creates a canonical dataset and starts its ingestion by invoking the step function
        If `size` is not provided, it will be inferred automatically
        """
        logger.info(
            {
                "message": "ingesting dataset",
                "collection_version_id": collection_version_id,
                "url": url,
                "existing_dataset_version_id": existing_dataset_version_id,
            }
        )
        if not self.uri_provider.validate(url):
            raise InvalidURIException(f"Trying to upload invalid URI: {url}")

        if file_size is None:
            file_info = self.uri_provider.get_file_info(url)
            file_size = file_info.size

        from backend.common.corpora_config import CorporaConfig

        max_file_size_gb = CorporaConfig().upload_max_file_size_gb * 2**30

        if file_size is not None and file_size > max_file_size_gb:
            raise MaxFileSizeExceededException(f"{url} exceeds the maximum allowed file size of {max_file_size_gb} Gb")

        # Ensure that the collection exists and is not published
        collection = self._assert_collection_version_unpublished(collection_version_id)

        # Creates a dataset version that the processing pipeline will point to
        new_dataset_version: DatasetVersion
        if existing_dataset_version_id is not None:
            # Ensure that the dataset belongs to the collection
            if existing_dataset_version_id not in [d.version_id for d in collection.datasets]:
                raise DatasetNotFoundException(
                    f"Dataset {existing_dataset_version_id.id} does not belong to the desired collection"
                )

            dataset_version = self.database_provider.get_dataset_version(existing_dataset_version_id)
            if dataset_version is None:
                raise DatasetNotFoundException(
                    f"Trying to replace non existent dataset {existing_dataset_version_id.id}"
                )

            if dataset_version.status.processing_status not in [
                DatasetProcessingStatus.SUCCESS,
                DatasetProcessingStatus.FAILURE,
                DatasetProcessingStatus.INITIALIZED,
            ]:
                raise DatasetInWrongStatusException(
                    f"Unable to reprocess dataset {existing_dataset_version_id.id}: processing status is "
                    f"{dataset_version.status.processing_status.name}"
                )

            # If a dataset is "empty", we will not replace it but instead reuse the existing version id.
            # Make sure that the conditions are not relaxed or this will break immutability.
            if (
                dataset_version.status.processing_status == DatasetProcessingStatus.INITIALIZED
                and not dataset_version.artifacts
                and dataset_version.canonical_dataset.published_at is None
            ):
                new_dataset_version = dataset_version
            else:
                new_dataset_version = self.database_provider.replace_dataset_in_collection_version(
                    collection_version_id, existing_dataset_version_id
                )
        else:
            new_dataset_version = self.database_provider.create_canonical_dataset(collection_version_id)
            # adds new dataset version to collection version
            self.database_provider.add_dataset_to_collection_version_mapping(
                collection_version_id, new_dataset_version.version_id
            )

        # Sets an initial processing status for the new dataset version
        self.database_provider.update_dataset_upload_status(new_dataset_version.version_id, DatasetUploadStatus.WAITING)
        self.database_provider.update_dataset_processing_status(
            new_dataset_version.version_id, DatasetProcessingStatus.INITIALIZED
        )

        # Starts the step function process
        if start_step_function:
            self.step_function_provider.start_step_function(collection_version_id, new_dataset_version.version_id, url)

        return (new_dataset_version.version_id, new_dataset_version.dataset_id)

    def remove_dataset_version(
        self,
        collection_version_id: CollectionVersionId,
        dataset_version_id: DatasetVersionId,
        delete_published: bool = False,
    ) -> None:
        """
        Removes a dataset version from an existing collection version
        """
        self._assert_collection_version_unpublished(collection_version_id)
        dataset_version = self.database_provider.get_dataset_version(dataset_version_id)
        if dataset_version.canonical_dataset.published_at and not delete_published:
            raise CollectionUpdateException from None
        if delete_published and not dataset_version.canonical_dataset.published_at:
            raise DatasetIsPrivateException from None
        self.database_provider.delete_dataset_from_collection_version(collection_version_id, dataset_version_id)

    def set_dataset_metadata(self, dataset_version_id: DatasetVersionId, metadata: DatasetMetadata) -> None:
        """
        Sets the metadata for a dataset version
        """
        self.database_provider.set_dataset_metadata(dataset_version_id, metadata)

    def get_all_mapped_datasets(self) -> List[DatasetVersion]:
        """
        Retrieves all the datasets from the database that belong to a published collection
        """
        datasets, _ = self.database_provider.get_all_mapped_datasets_and_collections()
        return datasets

    def get_datasets_for_collections(self, collections: Iterable[CollectionVersion]) -> Iterable[DatasetVersion]:
        datasets = []
        for collection in collections:
            datasest_ids = [d.id for d in collection.datasets]
            collection_datasets: List[DatasetVersion] = [
                self.database_provider.get_dataset_version(DatasetVersionId(id)) for id in datasest_ids
            ]
            datasets.extend(collection_datasets)
        return datasets

    def get_all_mapped_collection_versions_with_datasets(self) -> List[CollectionVersionWithPublishedDatasets]:
        """
        Retrieves all the datasets from the database that belong to a published collection
        """
        (
            mapped_dataset_versions,
            mapped_collection_versions,
        ) = self.database_provider.get_all_mapped_datasets_and_collections()
        return self._map_collection_version_to_published_dataset_versions(
            mapped_collection_versions, mapped_dataset_versions
        )

    def get_dataset_artifacts(self, dataset_version_id: DatasetVersionId) -> Iterable[DatasetArtifact]:
        """
        Returns all the artifacts for a dataset
        """
        dataset = self.database_provider.get_dataset_version(dataset_version_id)
        return dataset.artifacts

    def get_dataset_artifact_download_data(
        self, dataset_version_id: DatasetVersionId, artifact_id: DatasetArtifactId
    ) -> DatasetArtifactDownloadData:
        """
        Returns download data for an artifact, including a presigned URL
        """
        artifacts = self.get_dataset_artifacts(dataset_version_id)
        artifact = next((a for a in artifacts if a.id == artifact_id), None)

        if not artifact:
            raise ArtifactNotFoundException(f"Artifact {artifact_id} not found in dataset {dataset_version_id}")

        file_name = artifact.uri.split("/")[-1]
        file_type = artifact.type
        file_size = self.s3_provider.get_file_size(artifact.uri)
        presigned_url = self.s3_provider.generate_presigned_url(artifact.uri)

        return DatasetArtifactDownloadData(file_name, file_type, file_size, presigned_url)

    def get_dataset_status(self, dataset_version_id: DatasetVersionId) -> DatasetStatus:
        """
        Returns the dataset status for a specific dataset version
        """
        return self.database_provider.get_dataset_version_status(dataset_version_id)

    def update_dataset_version_status(
        self,
        dataset_version_id: DatasetVersionId,
        status_key: DatasetStatusKey,
        new_dataset_status: DatasetStatusGeneric,
        validation_message: Optional[str] = None,
    ) -> None:
        """
        TODO: split into two method, one for updating validation_message, and the other statuses.
        Updates the status of a dataset version.
        status_key can be one of: [upload, validation, cxg, rds, h5ad, processing]
        """
        if status_key == DatasetStatusKey.UPLOAD and isinstance(new_dataset_status, DatasetUploadStatus):
            self.database_provider.update_dataset_upload_status(dataset_version_id, new_dataset_status)
        elif status_key == DatasetStatusKey.PROCESSING and isinstance(new_dataset_status, DatasetProcessingStatus):
            self.database_provider.update_dataset_processing_status(dataset_version_id, new_dataset_status)
        elif status_key == DatasetStatusKey.VALIDATION and isinstance(new_dataset_status, DatasetValidationStatus):
            self.database_provider.update_dataset_validation_status(dataset_version_id, new_dataset_status)
        elif status_key == DatasetStatusKey.CXG and isinstance(new_dataset_status, DatasetConversionStatus):
            self.database_provider.update_dataset_conversion_status(
                dataset_version_id, "cxg_status", new_dataset_status
            )
        elif status_key == DatasetStatusKey.RDS and isinstance(new_dataset_status, DatasetConversionStatus):
            self.database_provider.update_dataset_conversion_status(
                dataset_version_id, "rds_status", new_dataset_status
            )
        elif status_key == DatasetStatusKey.H5AD and isinstance(new_dataset_status, DatasetConversionStatus):
            self.database_provider.update_dataset_conversion_status(
                dataset_version_id, "h5ad_status", new_dataset_status
            )
        else:
            raise DatasetUpdateException(
                f"Invalid status update for dataset {dataset_version_id}: cannot set {status_key} to "
                f"{new_dataset_status}"
            )

        if validation_message is not None:
            self.database_provider.update_dataset_validation_message(dataset_version_id, validation_message)

    def add_dataset_artifact(
        self, dataset_version_id: DatasetVersionId, artifact_type: str, artifact_uri: str
    ) -> DatasetArtifactId:
        """
        Registers an artifact to a dataset version.
        """

        # TODO: we should probably validate that artifact_uri is a valid S3 URI

        if artifact_type not in [artifact.value for artifact in DatasetArtifactType]:
            raise DatasetIngestException(f"Wrong artifact type for {dataset_version_id}: {artifact_type}")

        return self.database_provider.add_dataset_artifact(dataset_version_id, artifact_type, artifact_uri)

    def update_dataset_artifact(self, artifact_id: DatasetArtifactId, artifact_uri: str) -> None:
        """
        Updates uri for an existing artifact_id
        """
        self.database_provider.update_dataset_artifact(artifact_id, artifact_uri)

    def create_collection_version(self, collection_id: CollectionId) -> CollectionVersionWithDatasets:
        """
        Creates a collection version for an existing canonical collection.
        Also ensures that the collection does not have any active, unpublished version
        """

        all_versions = self.database_provider.get_all_versions_for_collection(collection_id)
        if not all_versions:
            raise CollectionVersionException(f"Collection {collection_id} cannot be found")

        if any(v for v in all_versions if v.published_at is None):
            raise CollectionVersionException(f"Collection {collection_id} already has an unpublished version")

        added_version_id = self.database_provider.add_collection_version(collection_id)
        return self.get_collection_version(added_version_id)

    def delete_collection_version(self, collection_version: CollectionVersionWithDatasets) -> None:
        """
        Deletes a collection version. This method will raise an error if the version is published.
        (Note: for performance reasons, the check is performed by the underlying layer)
        """
        unpublished_versions_of_published_datasets: List[DatasetVersion] = []
        unpublished_datasets: List[DatasetVersion] = []
        versions_to_delete_from_s3: List[DatasetVersion] = []
        for dv in collection_version.datasets:
            unpublished_versions = self.get_unpublished_dataset_versions(dv.canonical_dataset.dataset_id)
            versions_to_delete_from_s3.extend(unpublished_versions)  # All unpublished s3 assets are to be deleted
            if dv.canonical_dataset.published_at:
                unpublished_versions_of_published_datasets.extend(unpublished_versions)
            else:
                unpublished_datasets.append(dv.canonical_dataset)
        self.delete_dataset_version_assets(versions_to_delete_from_s3)

        # Delete DatasetVersionTable rows and DatasetArtifactTable rows if Dataset is published
        self.database_provider.delete_dataset_versions(unpublished_versions_of_published_datasets)

        # Delete DatasetTable, DatasetVersionTable, and DatasetArtifactTable rows if Dataset is unpublished (new)
        self.database_provider.delete_datasets(unpublished_datasets)

        # Delete the Collection version itself
        self.database_provider.delete_collection_version(collection_version.version_id)
        if not collection_version.canonical_collection.originally_published_at:
            # Collection was never published; delete CollectionTable row
            self.database_provider.delete_unpublished_collection(collection_version.collection_id)

    def delete_dataset_versions_from_public_bucket(self, dataset_version_ids: List[str]) -> List[str]:
        rdev_prefix = os.environ.get("REMOTE_DEV_PREFIX", "").strip("/")
        object_keys = set()
        for d_v_id in dataset_version_ids:
            for file_type in ("h5ad", "rds"):
                dataset_version_s3_object_key = f"{d_v_id}.{file_type}"
                if rdev_prefix:
                    dataset_version_s3_object_key = f"{rdev_prefix}/{dataset_version_s3_object_key}"
                object_keys.add(dataset_version_s3_object_key)
        self._delete_from_bucket(os.getenv("DATASETS_BUCKET"), list(object_keys))
        return list(object_keys)

    def delete_all_dataset_versions_from_public_bucket_for_collection(self, collection_id: CollectionId) -> List[str]:
        """
        Delete all associated publicly-accessible Datasets in s3
        """
        dataset_versions = self.database_provider.get_all_dataset_versions_for_collection(collection_id)
        return self.delete_dataset_versions_from_public_bucket([dv.version_id.id for dv in dataset_versions])

    def get_unpublished_dataset_versions(self, dataset_id: DatasetId) -> List[DatasetVersion]:
        """
        Get all unpublished versions for a Dataset that are currently in the database
        """
        all_versions = self.database_provider.get_all_versions_for_dataset(dataset_id)
        dataset = all_versions[0].canonical_dataset
        # Determine when the current mapped version was created
        mapped_version_created_at = (
            next(d_v.created_at for d_v in all_versions if d_v.version_id == dataset.dataset_version_id)
            if dataset.published_at
            else datetime.min
        )
        # Gather all Dataset versions that were created after the current mapped version
        return list(
            filter(
                lambda d_v: d_v.created_at > mapped_version_created_at,
                all_versions,
            )
        )

    def delete_dataset_version_assets(self, dataset_versions: List[DatasetVersion]) -> None:
        self.delete_dataset_versions_from_public_bucket([dv.version_id.id for dv in dataset_versions])
        self.delete_artifacts(reduce(lambda artifacts, dv: artifacts + dv.artifacts, dataset_versions, []))

    def tombstone_collection(self, collection_id: CollectionId) -> None:
        """
        Tombstone a published Collection. For admin use only.
        """
        self.delete_all_dataset_versions_from_public_bucket_for_collection(collection_id)
        self.database_provider.tombstone_collection(collection_id)

    def resurrect_collection(self, collection_id: CollectionId) -> None:
        """
        Resurrect a tombstoned Collection (untombstone Collection and un-delete public s3 assets). Doing so restores
        accessibility of all previous Collection versions and constituent Datasets **EXCEPT** for Datasets which were
        tombstoned individually before the Collection was tombstoned as a whole; such Datasets remain tombstoned after
        resurrection of the Collection. Effectively, this restores the Collection to its most recent state immediately
        prior to being tombstoned. For admin use only.
        """
        collection = self.get_canonical_collection(collection_id)
        if not collection.tombstoned:
            raise CollectionIsPublicException()

        # Individually-tombstoned Datasets must remain tombstoned after resurrecting Collection
        collection_versions = sorted(
            self.get_all_published_collection_versions_from_canonical(collection_id, get_tombstoned=True),
            key=lambda cv: cv.created_at,
        )
        dataset_ids_to_version_ids: Dict[str, List[str]] = defaultdict(list)
        tombstoned_datasets: Set[str] = set()  # Track individually-tombstoned Datasets
        previous_dataset_ids: Set[str] = set()
        for cv in collection_versions:
            [dataset_ids_to_version_ids[dv.dataset_id.id].append(dv.version_id.id) for dv in cv.datasets]
            current_dataset_ids: Set[str] = {dv.dataset_id.id for dv in cv.datasets}
            tombstoned_datasets.update(previous_dataset_ids.difference(current_dataset_ids))
            previous_dataset_ids = current_dataset_ids

        dataset_versions_to_restore: List[str] = []
        datasets_to_resurrect: List[str] = []
        for dataset_id, version_ids in dataset_ids_to_version_ids.items():
            if dataset_id not in tombstoned_datasets:
                dataset_versions_to_restore.extend(version_ids)
                datasets_to_resurrect.append(dataset_id)

        # Restore s3 public assets
        for dv_id in dataset_versions_to_restore:
            for ext in (DatasetArtifactType.H5AD, DatasetArtifactType.RDS):
                object_key = f"{dv_id}.{ext}"
                self.s3_provider.restore_object(os.getenv("DATASETS_BUCKET"), object_key)

        # Reset tombstone values in database
        self.database_provider.resurrect_collection(collection_id, datasets_to_resurrect)

    def publish_collection_version(self, version_id: CollectionVersionId) -> None:
        """
        Publishes a collection version.
        """
        version = self.database_provider.get_collection_version_with_datasets(version_id)

        if version.published_at is not None:
            raise CollectionPublishException("Cannot publish an already published collection")

        if len(version.datasets) == 0:
            raise CollectionPublishException("Cannot publish a collection with no datasets")

        schema_versions = {dataset.metadata.schema_version for dataset in version.datasets}
        if len(schema_versions) != 1:
            raise CollectionPublishException("Cannot publish a collection with datasets of different schema versions")
        schema_version = next(iter(schema_versions))

        date_of_last_publish = datetime.min
        has_dataset_revisions = False
        # if collection is a revision and has no changes to previous version's datasets--don't update 'revised_at'
        # used for cases where revision only contains collection-level metadata changes
        if version.canonical_collection.version_id is not None:
            date_of_last_publish = (
                version.canonical_collection.revised_at or version.canonical_collection.originally_published_at
            )
            canonical_version = self.database_provider.get_collection_version(version.canonical_collection.version_id)
            canonical_datasets = {dataset_version_id.id for dataset_version_id in canonical_version.datasets}
            version_datasets = {dataset.version_id.id for dataset in version.datasets}
            if canonical_datasets != version_datasets:
                has_dataset_revisions = True

        # Finalize Collection publication and delete any tombstoned assets
        dataset_version_ids_to_delete_from_s3 = self.database_provider.finalize_collection_version(
            version.collection_id, version_id, schema_version, update_revised_at=has_dataset_revisions
        )
        self.delete_dataset_versions_from_public_bucket(dataset_version_ids_to_delete_from_s3)

        # Handle cleanup of unpublished versions
        dataset_versions = self.database_provider.get_all_dataset_versions_for_collection(
            version.collection_id, from_date=date_of_last_publish
        )
        versions_to_delete = list(
            filter(lambda dv: dv.version_id.id not in {dv.version_id.id for dv in version.datasets}, dataset_versions)
        )
        self.delete_dataset_version_assets(versions_to_delete)
        self.database_provider.delete_dataset_versions(versions_to_delete)

    def get_dataset_version(self, dataset_version_id: DatasetVersionId) -> Optional[DatasetVersion]:
        """
        Returns a dataset version by id
        """
        return self.database_provider.get_dataset_version(dataset_version_id)

    def get_dataset_version_from_canonical(
        self, dataset_id: DatasetId, get_tombstoned: bool = False
    ) -> Optional[DatasetVersion]:
        """
        Given a canonical dataset id, returns its mapped dataset version, i.e. the dataset version
        that belongs to the most recently published collection. If never published, will return the most recent
        unpublished version.
        """
        dataset_version = self.database_provider.get_dataset_mapped_version(dataset_id, get_tombstoned)
        if dataset_version is not None:
            if dataset_version.canonical_dataset.tombstoned:
                raise DatasetIsTombstonedException()
            return dataset_version
        # Dataset has never been published
        return self.database_provider.get_most_recent_active_dataset_version(dataset_id)

    def get_prior_published_versions_for_dataset(self, dataset_id: DatasetId) -> List[PublishedDatasetVersion]:
        """
        Given a canonical dataset id, return all its DatasetVersions that have been part of published CollectionVersions
        """
        dataset_version = self.database_provider.get_dataset_mapped_version(dataset_id, get_tombstoned=True)
        if not dataset_version:
            return []
        if dataset_version.canonical_dataset.tombstoned:
            raise DatasetIsTombstonedException()
        collection_versions = self.database_provider.get_all_versions_for_collection(dataset_version.collection_id)
        if collection_versions[0].canonical_collection.tombstoned:
            raise DatasetIsTombstonedException()
        published_version_history = []
        found_version_ids = set()
        # sort to ensure we always find earliest instance of a dataset version first when iterating
        # needs None check as part of sort key to avoid TypeError on sorting list of datetimes + NoneTypes
        collection_versions = sorted(collection_versions, key=lambda cv: (cv.published_at is None, cv.published_at))
        for collection_version in collection_versions:
            # skip unpublished collection versions
            if collection_version.published_at is None:
                # per sorting, the remaining collection versions are unpublished
                break
            for dataset_version in collection_version.datasets:
                if (
                    dataset_version.dataset_id.id == dataset_id.id
                    and dataset_version.version_id.id not in found_version_ids
                ):
                    published_version = PublishedDatasetVersion(
                        collection_version_id=collection_version.version_id,
                        published_at=collection_version.published_at,
                        **vars(dataset_version),
                    )
                    found_version_ids.add(dataset_version.version_id.id)
                    published_version_history.append(published_version)
        return published_version_history

    def get_prior_published_dataset_version(self, dataset_version_id: DatasetVersionId) -> PublishedDatasetVersion:
        """
        Given a dataset version id, return the DatasetVersion, if it's been part of a published CollectionVersion
        """
        dataset_version = self.database_provider.get_dataset_version(dataset_version_id, get_tombstoned=True)
        if not dataset_version:
            return None
        if dataset_version.canonical_dataset.tombstoned:
            raise DatasetIsTombstonedException()
        collection_versions = self.database_provider.get_all_versions_for_collection(dataset_version.collection_id)
        try:
            published_at, collection_version_id = get_published_at_and_collection_version_id_else_not_found(
                dataset_version, collection_versions
            )
            return PublishedDatasetVersion(
                collection_version_id=collection_version_id,
                published_at=published_at,
                **vars(dataset_version),
            )
        except DatasetVersionNotFoundException:
            return None

    def _delete_from_bucket(self, bucket: str, keys: List[str] = None, prefix: str = None) -> None:
        try:
            if keys:
                self.s3_provider.delete_files(bucket, keys)
            if prefix:
                self.s3_provider.delete_prefix(bucket, prefix)
        except S3DeleteException as e:
            raise CollectionDeleteException("Attempt to delete public Datasets failed") from e

    def delete_artifacts(self, artifacts: List[DatasetArtifact]) -> None:
        for artifact in artifacts:
            matches_dict = S3_URI_REGEX.match(artifact.uri).groupdict()
            bucket, key, prefix = matches_dict["bucket"], matches_dict["key"], matches_dict["prefix"]
            self._delete_from_bucket(bucket, keys=[key] if key else None, prefix=prefix)

    def _get_collection_and_dataset(
        self, collection_id: str, dataset_id: str
    ) -> Tuple[CollectionVersionWithDatasets, DatasetVersion]:
        """
        Get collection and dataset by their ids. Will look up by both version and canonical id for both.
        """

        collection_version = self.get_collection_version_from_canonical(CollectionId(collection_id))
        if collection_version is None:
            collection_version = self.get_collection_version(CollectionVersionId(collection_id))
        if collection_version is None:
            raise CollectionNotFoundException()

        # Extract the dataset from the dataset list.
        try:
            dataset_version = next(
                d for d in collection_version.datasets if d.version_id.id == dataset_id or d.dataset_id.id == dataset_id
            )
        except StopIteration:
            raise DatasetNotFoundException() from None

        return collection_version, dataset_version

    def _map_collection_version_to_published_dataset_versions(
        self, mapped_collection_versions: List[CollectionVersion], mapped_dataset_versions: List[DatasetVersion]
    ) -> List[CollectionVersionWithPublishedDatasets]:
        """
        Given a list of collection versions and dataset versions, produce a list joining them into
        CollectionVersionWithPublishedDatasets objects
        """
        # Construct dict of collection_id: [Datasets]
        datasets_by_collection_id = defaultdict(list)
        [datasets_by_collection_id[d.collection_id.id].append(d) for d in mapped_dataset_versions]

        # Construct dict of collection_id: [all published Collection versions]
        all_collections_versions = self.database_provider.get_all_collections_versions()
        collection_versions_by_collection_id = defaultdict(list)
        [collection_versions_by_collection_id[c.collection_id.id].append(c) for c in all_collections_versions]

        # Determine published_at and collection_version_id for all mapped Dataset versions. Both values come from the
        # Collection version with which the Dataset version was first published.
        collections_with_published_datasets: List[CollectionVersionWithPublishedDatasets] = []
        for collection in mapped_collection_versions:
            mapped_dataset_versions_for_collection = datasets_by_collection_id[collection.collection_id.id]
            all_versions_for_collection = collection_versions_by_collection_id[collection.collection_id.id]
            published_datasets_for_collection: List[PublishedDatasetVersion] = []
            for mapped_dataset_version in mapped_dataset_versions_for_collection:
                published_at, collection_version_id = get_published_at_and_collection_version_id_else_not_found(
                    mapped_dataset_version, all_versions_for_collection
                )
                revised_at = (
                    None if mapped_dataset_version.canonical_dataset.published_at == published_at else published_at
                )
                mapped_dataset_version.canonical_dataset.revised_at = revised_at
                published_datasets_for_collection.append(
                    PublishedDatasetVersion(
                        collection_version_id=collection_version_id,
                        revised_at=revised_at,  # Duped logic to avoid pitfall; no effect in present implementation
                        published_at=mapped_dataset_version.canonical_dataset.published_at,
                        **vars(mapped_dataset_version),
                    )
                )
            collection.datasets = published_datasets_for_collection  # hack to allow unpacking via **vars() below
            collections_with_published_datasets.append(CollectionVersionWithPublishedDatasets(**vars(collection)))

        return collections_with_published_datasets

    def restore_previous_dataset_version(
        self, collection_version_id: CollectionVersionId, dataset_id: DatasetId
    ) -> None:
        """
        Restore the previous dataset version for a dataset.
        :param collection_version_id: The collection version to restore the dataset version. It must be in a mutable
        state
        :param dataset_id: The dataset id to restore the previous version of.
        """

        collection = self.database_provider.get_collection_version_with_datasets(collection_version_id)
        if collection.is_published():
            raise CollectionIsPublishedException(
                f"Collection {collection_version_id} is published, unable to restore previous dataset version"
            )
        current_version: DatasetVersion = [dv for dv in collection.datasets if dv.dataset_id == dataset_id][0]
        previous_version_id: DatasetVersionId = self.database_provider.get_previous_dataset_version_id(dataset_id)
        if previous_version_id is None:
            raise NoPreviousDatasetVersionException(f"No previous dataset version for dataset {dataset_id.id}")
        logger.info(
            {
                "message": "Restoring previous dataset version",
                "dataset_id": dataset_id.id,
                "replace_version_id": current_version.version_id.id,
                "restored_version_id": previous_version_id.id,
            }
        )
        self.database_provider.replace_dataset_in_collection_version(
            collection_version_id, current_version.version_id, previous_version_id
        )
