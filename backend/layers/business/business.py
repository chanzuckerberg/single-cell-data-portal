from dataclasses import dataclass
from typing import Iterable, Optional
from backend.corpora.common.providers.crossref_provider import CrossrefDOINotFoundException, CrossrefException
from backend.layers.business.entities import CollectionMetadataUpdate, CollectionQueryFilter, DatasetArtifactDownloadData
from backend.layers.business.exceptions import CollectionCreationException, CollectionUpdateException, DatasetIngestException

from backend.layers.common.entities import (
    CollectionId,
    CollectionLinkType,
    CollectionMetadata,
    CollectionVersion,
    CollectionVersionId,
    DatasetArtifact,
    DatasetId,
    DatasetMetadata,
    DatasetProcessingStatus,
    DatasetStatus,
    DatasetStatusGeneric,
    DatasetVersion,
    DatasetVersionId,
    Link,
)
from typing import Iterable, List, Optional

from backend.layers.common.entities import CollectionId, CollectionMetadata, CollectionVersion, CollectionVersionId, DatasetArtifact, DatasetId, DatasetStatus, DatasetVersion, DatasetVersionId, Link
from backend.layers.persistence.persistence import DatabaseProviderInterface
from backend.layers.thirdparty.crossref_provider import CrossrefProviderInterface
from backend.layers.thirdparty.step_function_provider import StepFunctionProviderInterface

import re
from urllib.parse import urlparse
from backend.layers.common import validation
import logging


class BusinessLogicInterface:

    # Get_collections
    # Replaces get_collections_list and get_collections_index
    # Accepts a CollectionFilter class (or kwargs) with:
    # Date bounds
    # Visibility
    # Ownership
    # List of fields to be returned
    # Returns a list of dictionaries that only includes the selected fields for each collection
    # It should NOT add any information that is required by the current API for compatibility reasons (e.g. access_write). These will be delegated to the upper layer
    # It should NOT do any operation that is required by Curation API assumptions (e.g. remove None values)

    def get_collections(self, filter: CollectionQueryFilter) -> Iterable[CollectionVersion]:
        pass

    # Get_collection
    # Replaces get_collection_details
    # Returns a single collection, with no filtering options
    # Accepts:
    # Collection_id
    # Should reuse most of of the code from the method above

    def get_published_collection_version(self, collection_id: CollectionId) -> CollectionVersion:
        pass

    def get_collection_version(self, version_id: CollectionVersionId) -> CollectionVersion:
        pass

    # Create_collection
    # Replaces the current create_collection
    # Accepts:
    # A dictionary (or Class) with the body of the collection to be created
    # Should validate the body accepted as param (see existing verify_collection_body)
    # Should call CrossrefProvider to retrieve publisher metadata information
    # This method currently collects errors in a list, which will be piped upstream to the API response. This is a good idea but it should be refactor into a generalized pattern (otherwise we’ll “pollute” the business layer with logic specific to the API layer).

    def create_collection(self, collection_metadata: CollectionMetadata) -> CollectionVersion:
        pass

    # Delete_collection
    # Replaces the current delete_collection
    # Accepts:
    # Collection_id
    # Performs authorization on user/collection

    def delete_collection(self, collection_id: CollectionId) -> None:
        pass

    # Update_collection
    # Replaces the current update_collection
    # Accepts:
    # Collection_id
    # A dataclass with the body to be updated
    # Should validate the body
    # Should handle DOI updates (re-use the existing logic with minimal refactors)
    # Can either return nothing or the metadata of the updated collection

    # TODO: body should be a dataclass?
    def update_collection_version(self, version_id: CollectionVersionId, body: CollectionMetadataUpdate) -> None:
        pass

    # Create_collection_version
    # Replaces the current post_collection_revision
    # Accepts:
    # Collection_id
    # Performs authorization on the collection
    # Since revision logic is database specific, it delegates to the underlying layer
    # Returns a handle to the revised collection (either id or the full collection metadata)

    def create_collection_version(self, collection_id: CollectionId) -> CollectionVersion:
        pass

    def delete_collection_version(self, version_id: CollectionVersionId) -> None:
        pass

    # Publish_collection
    # Replaces post (in publish.py)
    # Accepts:
    # Collection_id
    # Performs validation to make sure that the collection can be published
    # [Currently] accepts data_submission_policy_version: what is this for?
    # [Currently] triggers Cloudfront invalidation for the index endpoints. This should arguably NOT be done here but by the API layer
    # Since revision logic is database specific, it delegates to the underlying layer

    def publish_collection_version(self, version_id: CollectionVersionId) -> None:
        pass

    # Ingest_dataset
    # Replaces the existing Upload_from_link
    # Potentially, also replaces relink (I am not sure why they are 2 separate functions)
    # Accepts:
    # Collection_id
    # URL of the uploadable dataset
    # [Optional] a dataset_id to be replaced
    # This is one of the most complex functions. Other than the database provider, It will need two additional providers:
    # StepFunctionProvider (to call the SFN that triggers the upload)
    # DropboxProvider to interface with Dropbox (could be more generic: RemoteFileProvider?)
    # Should handle exceptions from all providers:
    # Should only raise custom exceptions

    def ingest_dataset(
        self,
        collection_version_id: CollectionVersionId,
        url: str,
        existing_dataset_version_id: Optional[DatasetVersionId],
    ) -> DatasetVersionId:
        pass

    # Get_all_datasets
    # Replaces get_dataset_index

    def get_all_datasets(self) -> Iterable[DatasetVersion]:
        pass

    # Delete_dataset
    # Replaces delete_dataset

    def delete_dataset(self, dataset_version_id: DatasetVersionId) -> None:
        pass

    # get_dataset_assets
    # Replaces get_dataset_assets

    def get_dataset_artifacts(self, dataset_id: DatasetId) -> Iterable[DatasetArtifact]:
        pass

    # Download_dataset_asset
    # Replaces post_dataset_asset

    def get_dataset_artifact_download_data(
        self, dataset_id: DatasetId, artifact_id: str
    ) -> DatasetArtifactDownloadData:
        pass

    def update_dataset_version_status(
        self, dataset_version_id: DatasetVersionId, new_dataset_status: DatasetStatusGeneric
    ) -> None:
        pass

    def add_dataset_artifact(self, dataset_version_id: DatasetVersionId, artifact_type: str, artifact_uri: str) -> None:
        pass

    # Get_dataset_status
    # Replaces get_status

    def get_dataset_status(self, dataset_id: DatasetId) -> DatasetStatus:
        pass


# TODO: move it to a separate file
class BusinessLogic(BusinessLogicInterface):

    database_provider: DatabaseProviderInterface
    crossref_provider: CrossrefProviderInterface
    step_function_provider: StepFunctionProviderInterface

    def __init__(
        self,
        database_provider: DatabaseProviderInterface,
        crossref_provider: CrossrefProviderInterface,
        step_function_provider: StepFunctionProviderInterface,
    ) -> None:
        self.crossref_provider = crossref_provider
        self.database_provider = database_provider
        self.step_function_provider = step_function_provider
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
        doi = next((link.uri for link in collection_metadata.links if link.type == "doi"), None)

        if doi is not None:
            publisher_metadata = self._get_publisher_metadata(doi, errors)
        else:
            publisher_metadata = None

        if errors:
            raise CollectionCreationException(errors)

        created_version = self.database_provider.create_canonical_collection(owner, collection_metadata)

        # TODO: can collapse with `create_canonical_collection`
        if publisher_metadata:
            self.database_provider.save_collection_publisher_metadata(created_version.version_id, publisher_metadata)

        return created_version

    def get_published_collection_version(self, collection_id: CollectionId) -> Optional[CollectionVersion]:
        """
        Returns the published collection version that belongs to a canonical collection.
        Returns None if no published collection exists
        """
        return self.database_provider.get_collection_mapped_version(collection_id)

    def get_collection_version(self, version_id: CollectionVersionId) -> CollectionVersion:
        """
        Returns a specific collection version
        """
        return self.database_provider.get_collection_version(version_id)

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
        old_doi = next((link.uri for link in current_version.metadata.links if link.type == "doi"), None)
        if body.links is None:
            new_doi = None
        else:
            new_doi = next((link.uri for link in body.links if link.type == "doi"), None)

        if old_doi and new_doi is None:
            # If the DOI was deleted, remove the publisher_metadata field
            self.database_provider.save_collection_publisher_metadata(version_id, None)
        elif (new_doi is not None) and new_doi != old_doi:
            # If the DOI has changed, fetch and update the metadata
            publisher_metadata = self._get_publisher_metadata(new_doi, errors)
            self.database_provider.save_collection_publisher_metadata(version_id, publisher_metadata)

        if errors:
            raise CollectionUpdateException(errors)

        # Merge the updated fields in the existing metadata object
        current_metadata = current_version.metadata
        for field in vars(body):
            if hasattr(body, field):
                setattr(current_metadata, field, getattr(body, field))

        self.database_provider.save_collection_metadata(version_id, current_metadata)

    def ingest_dataset(
        self,
        collection_version_id: CollectionVersionId,
        url: str,
        existing_dataset_version_id: Optional[DatasetVersionId]) -> DatasetVersionId:
        """
        Creates a canonical dataset and starts its ingestion by invoking the step function
        """

        # TODO: add link validation

        # Verify Dropbox URL
        # valid_link = from_url(url)
        # if not valid_link:
        #     raise InvalidParametersHTTPException(detail="The dropbox shared link is invalid.")

        # TODO: add file_info through a provider
        # Get file info
        # try:
        #     resp = valid_link.file_info()
        # except requests.HTTPError:
        #     raise InvalidParametersHTTPException(detail="The URL provided causes an error with Dropbox.")
        # except MissingHeaderException as ex:
        #     raise InvalidParametersHTTPException(detail=ex.detail)

        # file_size = resp.get("size")

        collection_version = self.database_provider.get_collection_version(collection_version_id)
        if collection_version is None:
            raise DatasetIngestException(f"Collection version {collection_version_id} does not exist")


        new_dataset_version: DatasetVersion

        if existing_dataset_version_id is not None:
            dataset_version = self.database_provider.get_dataset_version(existing_dataset_version_id)
            if dataset_version is None:
                raise DatasetIngestException(f"Trying to replace non existant dataset {existing_dataset_version_id}")

            if dataset_version.status.processing_status not in [
                DatasetProcessingStatus.SUCCESS,
                DatasetProcessingStatus.FAILURE,
                DatasetProcessingStatus.INITIALIZED,
            ]:
                raise DatasetIngestException(
                    f"Unable to reprocess dataset {existing_dataset_version_id}: {dataset_version.status.processing_status=}"
                )

            # TODO: `add_dataset_version` should not take metadata, since it will be replaced by the processing pipeline
            new_dataset_version = self.database_provider.add_dataset_version(dataset_version.dataset_id, dataset_version.metadata)
        else:
            self.database_provider.create_canonical_dataset(collection_version.version_id)

        # if dataset:
        #     # Update dataset
        #     if dataset.processing_status.processing_status not in [
        #         ProcessingStatus.SUCCESS,
        #         ProcessingStatus.FAILURE,
        #         ProcessingStatus.INITIALIZED,
        #     ]:
        #         raise InvalidProcessingStateException(
        #             f"Unable to reprocess dataset {dataset_id}: {dataset.processing_status.processing_status=}"
        #         )
        #     else:
        #         dataset.reprocess()

        # else:
        #     # Add new dataset
        #     dataset = Dataset.create(db_session, collection=collection)

        # dataset.update(processing_status=dataset.new_processing_status())

        # # Start processing link
        # start_upload_sfn(collection_id, dataset.id, url)

        # return dataset.id


    def get_all_datasets(self) -> Iterable[DatasetVersion]:
        pass

    def delete_dataset(self, dataset_version_id: DatasetVersionId) -> None:
        pass

    def set_dataset_metadata(self, dataset_version_id: DatasetVersionId, metadata: DatasetMetadata) -> None:
        pass




    # def start_upload_sfn(collection_id, dataset_id, url):
    #     input_parameters = {
    #         "collection_id": collection_id,
    #         "url": url,
    #         "dataset_id": dataset_id,
    #     }
    #     sfn_name = f"{dataset_id}_{int(time.time())}"
    #     response = get_stepfunctions_client().start_execution(
    #         stateMachineArn=CorporaConfig().upload_sfn_arn,
    #         name=sfn_name,
    #         input=json.dumps(input_parameters),
    #     )
    #     return response



    # def upload_from_link(collection_id: str, token_info: dict, url: str, dataset_id: str = None):
    #     db_session = g.db_session

    #     # Verify Dropbox URL
    #     valid_link = from_url(url)
    #     if not valid_link:
    #         raise InvalidParametersHTTPException(detail="The dropbox shared link is invalid.")

    #     # Get file info
    #     try:
    #         resp = valid_link.file_info()
    #     except requests.HTTPError:
    #         raise InvalidParametersHTTPException(detail="The URL provided causes an error with Dropbox.")
    #     except MissingHeaderException as ex:
    #         raise InvalidParametersHTTPException(detail=ex.detail)

    #     file_size = resp.get("size")

    #     try:
    #         return upload(
    #             db_session,
    #             collection_id=collection_id,
    #             url=url,
    #             file_size=file_size,
    #             user=token_info["sub"],
    #             scope=token_info["scope"],
    #             dataset_id=dataset_id,
    #         )
    #     except MaxFileSizeExceededException:
    #         raise TooLargeHTTPException()
    #     except InvalidFileFormatException:
    #         raise InvalidParametersHTTPException(detail="The file referred to by the link is not a support file format.")
    #     except NonExistentCollectionException:
    #         raise ForbiddenHTTPException()
    #     except InvalidProcessingStateException:
    #         raise MethodNotAllowedException(
    #             detail="Submission failed. A dataset cannot be updated while a previous update for the same dataset is in "
    #             "progress. Please cancel the current submission by deleting the dataset, or wait until the submission has "
    #             "finished processing.",
    #         )
    #     except NonExistentDatasetException:
    #         raise NotFoundHTTPException()


    # def upload(
    #     db_session: Session,
    #     collection_id: str,
    #     url: str,
    #     file_size: int,
    #     user: str,
    #     scope: str = None,
    #     dataset_id: str = None,
    # ) -> str:
    #     max_file_size_gb = CorporaConfig().upload_max_file_size_gb * GB
    #     if file_size is not None and file_size > max_file_size_gb:
    #         raise MaxFileSizeExceededException(f"{url} exceeds the maximum allowed file size of {max_file_size_gb} Gb")

    #     # Check if datasets can be added to the collection
    #     collection = Collection.get_collection(
    #         db_session,
    #         collection_id,
    #         visibility=CollectionVisibility.PRIVATE,  # Do not allow changes to public Collections
    #         owner=owner_or_allowed(user, scope) if scope else user,
    #     )
    #     if not collection:
    #         raise NonExistentCollectionException(f"Collection {collection_id} does not exist")

    #     # Check if a dataset already exists
    #     if dataset_id:
    #         dataset = Dataset.get(db_session, dataset_id, collection_id=collection_id)
    #         if not dataset:
    #             raise NonExistentDatasetException(f"Dataset {dataset_id} does not exist")
    #     else:
    #         dataset = None

    #     if dataset:
    #         # Update dataset
    #         if dataset.processing_status.processing_status not in [
    #             ProcessingStatus.SUCCESS,
    #             ProcessingStatus.FAILURE,
    #             ProcessingStatus.INITIALIZED,
    #         ]:
    #             raise InvalidProcessingStateException(
    #                 f"Unable to reprocess dataset {dataset_id}: {dataset.processing_status.processing_status=}"
    #             )
    #         else:
    #             dataset.reprocess()

    #     else:
    #         # Add new dataset
    #         dataset = Dataset.create(db_session, collection=collection)

    #     dataset.update(processing_status=dataset.new_processing_status())

    #     # Start processing link
    #     start_upload_sfn(collection_id, dataset.id, url)

    #     return dataset.id
