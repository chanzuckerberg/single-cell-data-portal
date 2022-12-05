from dataclasses import dataclass
from typing import Iterable, Optional, Tuple
from backend.common.providers.crossref_provider import CrossrefDOINotFoundException, CrossrefException
from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.business.entities import CollectionMetadataUpdate, CollectionQueryFilter, DatasetArtifactDownloadData
from backend.layers.business.exceptions import ArtifactNotFoundException, CollectionCreationException, CollectionIsPublishedException, CollectionNotFoundException, CollectionPublishException, CollectionUpdateException, CollectionVersionException, DatasetInWrongStatusException, DatasetIngestException, DatasetNotFoundException, DatasetUpdateException, InvalidURIException, MaxFileSizeExceededException
from backend.layers.common.entities import (
    CollectionId,
    CollectionLinkType,
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
    DatasetStatusGeneric,
    DatasetStatusKey,
    DatasetUploadStatus,
    DatasetValidationStatus,
    DatasetVersion,
    DatasetVersionId,
    Link,
)
from typing import Iterable, List, Optional

import copy

from backend.layers.persistence.persistence_interface import DatabaseProviderInterface
from backend.layers.thirdparty.crossref_provider import CrossrefProviderInterface
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.step_function_provider import StepFunctionProviderInterface

import re
from urllib.parse import urlparse
from backend.layers.common import validation
import logging

from backend.layers.thirdparty.uri_provider import UriProvider, UriProviderInterface


class BusinessLogic(BusinessLogicInterface):

    database_provider: DatabaseProviderInterface
    crossref_provider: CrossrefProviderInterface
    step_function_provider: StepFunctionProviderInterface
    s3_provider: S3Provider
    uri_provider: UriProviderInterface

    def __init__(
        self,
        database_provider: DatabaseProviderInterface,
        crossref_provider: CrossrefProviderInterface,
        step_function_provider: StepFunctionProviderInterface,
        s3_provider: S3Provider,
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

    def create_collection(self, owner: str, collection_metadata: CollectionMetadata) -> CollectionVersion:
        """
        Creates a collection using the specified metadata. If a DOI is defined, will also
        retrieve publisher metadata from Crossref and add it to the collection.
        """

        errors = []
        validation.verify_collection_metadata(collection_metadata, errors)

        # TODO: Maybe switch link.type to be an enum
        doi = next((link.uri for link in collection_metadata.links if link.type == "DOI"), None)
        print("doi", doi)


        if doi is not None:
            publisher_metadata = self._get_publisher_metadata(doi, errors)
        else:
            publisher_metadata = None

        print(errors)

        if errors:
            raise CollectionCreationException(errors)

        created_version = self.database_provider.create_canonical_collection(owner, collection_metadata)

        # TODO: can collapse with `create_canonical_collection`
        if publisher_metadata:
            self.database_provider.save_collection_publisher_metadata(created_version.version_id, publisher_metadata)
            # Add this to the returned object, to save another `get_collection` call
            created_version.publisher_metadata = publisher_metadata 

        return created_version

    def get_published_collection_version(self, collection_id: CollectionId) -> Optional[CollectionVersion]:
        """
        Returns the published collection version that belongs to a canonical collection.
        Returns None if no published collection exists
        """
        return self.database_provider.get_collection_mapped_version(collection_id)

    def get_collection_version(self, version_id: CollectionVersionId) -> CollectionVersionWithDatasets:
        """
        Returns a specific collection version
        """
        return self.database_provider.get_collection_version_with_datasets(version_id)

    def get_collection_versions_from_canonical(self, collection_id: CollectionId) -> Iterable[CollectionVersion]:
        """
        Returns all the collection versions connected to a canonical collection
        """
        return self.database_provider.get_all_versions_for_collection(collection_id)

    def get_collection_version_from_canonical(self, collection_id: CollectionId) -> Optional[CollectionVersion]:
        """
        Returns the published collection version mapped to a canonical collection, if available.
        Otherwise will return the active unpublished version.
        """
        published_version = self.get_published_collection_version(collection_id)
        if published_version is not None:
            return published_version
        else:
            all_versions = self.get_collection_versions_from_canonical(collection_id)
            return next(v for v in all_versions if v.published_at is None)

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
                filter.is_published is None or 
                (filter.is_published is True and version.published_at is not None) or 
                (filter.is_published is False and version.published_at is None)
            )
            owner = (
                filter.owner is None or filter.owner == version.owner
            )
            return published and owner

        for collection_version in iterable:
            if predicate(collection_version):
                yield collection_version

    def update_collection_version(self, version_id: CollectionVersionId, body: CollectionMetadataUpdate) -> None:
        """
        Updates a collection version by replacing parts of its metadata. If the DOI in the links changed,
        it should also update its publisher metadata
        """

        # TODO: we could collect all the changes and issue a single update at the end
        # TODO: CollectionMetadataUpdate should probably be used for collection creation as well
        # TODO: link.type should DEFINITELY move to an enum. pylance will help with the refactor

        errors = []
        validation.verify_collection_metadata_update(body, errors)

        current_version = self.get_collection_version(version_id)
        if current_version.published_at is not None:
            raise CollectionUpdateException(["Cannot update a published collection"])
        
        # Determine if the DOI has changed
        old_doi = next((link.uri for link in current_version.metadata.links if link.type == "DOI"), None)
        if body.links is None:
            new_doi = None
        else:
            new_doi = next((link.uri for link in body.links if link.type == "DOI"), None)

        if old_doi and new_doi is None:
            # If the DOI was deleted, remove the publisher_metadata field
            self.database_provider.save_collection_publisher_metadata(version_id, None)
        elif (new_doi is not None) and new_doi != old_doi:
            # If the DOI has changed, fetch and update the metadata
            publisher_metadata = self._get_publisher_metadata(new_doi, errors)
            self.database_provider.save_collection_publisher_metadata(version_id, publisher_metadata)

        if errors:
            raise CollectionUpdateException(errors)

        # Merge the updated fields in the existing metadata object. Use a copy to ensure immutability.
        new_metadata = copy.deepcopy(current_version.metadata)
        for field in vars(body):
            if hasattr(body, field) and (value := getattr(body, field)) is not None:
                setattr(new_metadata, field, value)

        self.database_provider.save_collection_metadata(version_id, new_metadata)

    def _assert_collection_version_unpublished(self, collection_version_id: CollectionVersionId) -> CollectionVersionWithDatasets:
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

    # TODO: Alternatives: 1) return DatasetVersion 2) Return a new class
    def ingest_dataset(
        self,
        collection_version_id: CollectionVersionId,
        url: str,
        existing_dataset_version_id: Optional[DatasetVersionId]) -> Tuple[DatasetVersionId, DatasetId]:
        """
        Creates a canonical dataset and starts its ingestion by invoking the step function
        """

        if not self.uri_provider.validate(url):
            raise InvalidURIException(f"Trying to upload invalid URI: {url}")

        file_info = self.uri_provider.get_file_info(url)

        max_file_size_gb = 30 * 2**30 # TODO: read it from the config - requires smart mocking
        # max_file_size_gb = CorporaConfig().upload_max_file_size_gb * GB

        if file_info.size is not None and file_info.size > max_file_size_gb:
            raise MaxFileSizeExceededException(f"{url} exceeds the maximum allowed file size of {max_file_size_gb} Gb")

        # Ensure that the collection exists and is not published
        collection = self._assert_collection_version_unpublished(collection_version_id)

        # Creates a dataset version that the processing pipeline will point to
        new_dataset_version: DatasetVersion
        if existing_dataset_version_id is not None:
            # Ensure that the dataset belongs to the collection
            if existing_dataset_version_id not in [d.version_id for d in collection.datasets]:
                raise DatasetNotFoundException(f"Dataset {existing_dataset_version_id.id} does not belong to the desired collection")

            dataset_version = self.database_provider.get_dataset_version(existing_dataset_version_id)
            if dataset_version is None:
                raise DatasetNotFoundException(f"Trying to replace non existent dataset {existing_dataset_version_id.id}")

            if dataset_version.status.processing_status not in [
                DatasetProcessingStatus.SUCCESS,
                DatasetProcessingStatus.FAILURE,
                DatasetProcessingStatus.INITIALIZED,
            ]:
                raise DatasetInWrongStatusException(
                    f"Unable to reprocess dataset {existing_dataset_version_id.id}: processing status is {dataset_version.status.processing_status.name}"
                )

            # TODO: this method could very well be called `add_dataset_version`
            new_dataset_version = self.database_provider.replace_dataset_in_collection_version(collection_version_id, existing_dataset_version_id)
        else:
            new_dataset_version = self.database_provider.create_canonical_dataset(collection_version_id)
            # adds new dataset version to collection version
            self.database_provider.add_dataset_to_collection_version_mapping(collection_version_id, new_dataset_version.version_id)

        # Sets an initial processing status for the new dataset version
        self.database_provider.update_dataset_upload_status(new_dataset_version.version_id, DatasetUploadStatus.WAITING)
        self.database_provider.update_dataset_processing_status(new_dataset_version.version_id, DatasetProcessingStatus.PENDING)

        # Starts the step function process
        self.step_function_provider.start_step_function(collection_version_id, new_dataset_version.version_id, url)

        return (new_dataset_version.version_id, new_dataset_version.dataset_id)


    def remove_dataset_version(self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId) -> None:
        """
        Removes a dataset version from an existing collection version
        """
        self._assert_collection_version_unpublished(collection_version_id)
        self.database_provider.delete_dataset_from_collection_version(collection_version_id, dataset_version_id)

    def set_dataset_metadata(self, dataset_version_id: DatasetVersionId, metadata: DatasetMetadata) -> None:
        """
        Sets the metadata for a dataset version
        """
        self.database_provider.set_dataset_metadata(dataset_version_id, metadata)

    def get_all_published_datasets(self) -> Iterable[DatasetVersion]:
        """
        Retrieves all the datasets from the database that belong to a published collection
        """
        return self.database_provider.get_all_datasets()

    def get_dataset_artifacts(self, dataset_version_id: DatasetVersionId) -> Iterable[DatasetArtifact]:
        """
        Returns all the artifacts for a dataset
        """
        dataset = self.database_provider.get_dataset_version(dataset_version_id)
        return dataset.artifacts

    def get_dataset_artifact_download_data(self, dataset_version_id: DatasetVersionId, artifact_id: DatasetArtifactId) -> DatasetArtifactDownloadData:
        """
        Returns download data for an artifact, including a presigned URL
        """
        artifacts = self.get_dataset_artifacts(dataset_version_id)
        artifact = next((a for a in artifacts if a.id == artifact_id), None)

        if not artifact:
            raise ArtifactNotFoundException(f"Artifact {artifact_id} not found in dataset {dataset_version_id}")

        artifact_url = urlparse(artifact.uri)

        file_name = artifact_url.path[1:]
        file_type = artifact.type
        file_size = self.s3_provider.get_file_size(artifact.uri)
        presigned_url = self.s3_provider.generate_presigned_url(artifact.uri)

        return DatasetArtifactDownloadData(file_name, file_type, file_size, presigned_url)

    def get_dataset_status(self, dataset_version_id: DatasetVersionId) -> DatasetStatus:
        """
        Returns the dataset status for a specific dataset version
        """
        return self.database_provider.get_dataset_version_status(dataset_version_id)

    def update_dataset_version_status(self, dataset_version_id: DatasetVersionId, status_key: DatasetStatusKey, new_dataset_status: DatasetStatusGeneric) -> None:
        """
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
            self.database_provider.update_dataset_conversion_status(dataset_version_id, "cxg_status", new_dataset_status)
        elif status_key == DatasetStatusKey.RDS and isinstance(new_dataset_status, DatasetConversionStatus):
            self.database_provider.update_dataset_conversion_status(dataset_version_id, "rds_status", new_dataset_status)
        elif status_key == DatasetStatusKey.H5AD and isinstance(new_dataset_status, DatasetConversionStatus):
            self.database_provider.update_dataset_conversion_status(dataset_version_id, "h5ad_status", new_dataset_status)
        else:
            raise DatasetUpdateException(f"Invalid status update for dataset {dataset_version_id}: cannot set {status_key} to {new_dataset_status}")

    def add_dataset_artifact(self, dataset_version_id: DatasetVersionId, artifact_type: str, artifact_uri: str) -> DatasetArtifactId:
        """
        Registers an artifact to a dataset version.
        """
        
        # TODO: we should probably validate that artifact_uri is a valid S3 URI

        if artifact_type not in ["H5AD", "CXG", "RDS"]:
            raise DatasetIngestException(f"Wrong artifact type for {dataset_version_id}: {artifact_type}")

        return self.database_provider.add_dataset_artifact(dataset_version_id, artifact_type, artifact_uri)

    def create_collection_version(self, collection_id: CollectionId) -> CollectionVersionWithDatasets:
        """
        Creates a collection version for an existing canonical collection.
        Also ensures that the collection does not have any active, unpublished version
        """

        try:
            all_versions = self.database_provider.get_all_versions_for_collection(collection_id)
        except Exception:  # TODO: maybe add a RecordNotFound exception for finer grained exceptions
            raise CollectionVersionException(f"Collection {collection_id} cannot be found")

        if any(v for v in all_versions if v.published_at is None):
            raise CollectionVersionException(f"Collection {collection_id} already has an unpublished version")

        added_version_id = self.database_provider.add_collection_version(collection_id)
        return self.get_collection_version(added_version_id)

    def delete_collection_version(self, version_id: CollectionVersionId) -> None:
        self.database_provider.delete_collection_version(version_id)

    def publish_collection_version(self, version_id: CollectionVersionId) -> None:
        version = self.database_provider.get_collection_version(version_id)

        if version.published_at is not None:
            raise CollectionPublishException("Cannot publish an already published collection")

        if len(version.datasets) == 0:
            raise CollectionPublishException("Cannot publish a collection with no datasets")

        self.database_provider.finalize_collection_version(version.collection_id, version_id)

    def get_dataset_version(self, dataset_version_id: DatasetVersionId) -> Optional[DatasetVersion]:
        return self.database_provider.get_dataset_version(dataset_version_id)

    def get_dataset_version_from_canonical(self, dataset_id: DatasetId) -> Optional[DatasetVersion]:
        return self.database_provider.get_dataset_mapped_version(dataset_id)
