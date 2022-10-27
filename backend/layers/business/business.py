from dataclasses import dataclass
from typing import Iterable, Optional

from backend.layers.common.entities import CollectionMetadata, CollectionVersion, DatasetArtifact, DatasetStatus, DatasetVersion
from backend.layers.persistence.persistence import DatabaseProviderInterface
from backend.layers.thirdparty.crossref_provider import CrossrefProviderInterface
from backend.layers.thirdparty.step_function_provider import StepFunctionProviderInterface


@dataclass
class CollectionQueryFilter:
    visibility: Optional[str] # Use an enum
    owner: Optional[str]
    # TODO: add list of fields to be returned (if needed)

@dataclass
class UserInfo:
    user_id: str
    token: str


class BusinessLogicInterface:

    # Get_collections
    # Replaces get_collections_list and get_collections_index
    # Accepts a CollectionFilter class (or kwargs) with:
    # Date bounds
    # Visibility
    # Ownership
    # List of fields to be returned
    # Accepts an UserInfo object
    # Returns a list of dictionaries that only includes the selected fields for each collection
    # It should NOT add any information that is required by the current API for compatibility reasons (e.g. access_write). These will be delegated to the upper layer
    # It should NOT do any operation that is required by Curation API assumptions (e.g. remove None values)

    def get_collections(self, filter: CollectionQueryFilter, user_info: UserInfo) -> Iterable[CollectionVersion]:
        pass


    # Get_collection
    # Replaces get_collection_details
    # Returns a single collection, with no filtering options
    # Accepts:
    # Collection_id
    # UserInfo for validation
    # Should reuse most of of the code from the method above

    def get_collection(self, collection_id: str, user_info: UserInfo) -> CollectionVersion:
        pass

    def get_collection_version(self, version_id: str, user_info: UserInfo) -> CollectionVersion
        pass

    # Create_collection
    # Replaces the current create_collection
    # Accepts:
    # A dictionary (or Class) with the body of the collection to be created
    # UserInfo
    # Should validate the body accepted as param (see existing verify_collection_body)
    # Should call CrossrefProvider to retrieve publisher metadata information
    # This method currently collects errors in a list, which will be piped upstream to the API response. This is a good idea but it should be refactor into a generalized pattern (otherwise we’ll “pollute” the business layer with logic specific to the API layer).

    def create_collection(self, collection_metadata: CollectionMetadata, user_info: UserInfo) -> CollectionVersion:
        pass

    # Delete_collection
    # Replaces the current delete_collection
    # Accepts:
    # Collection_id
    # UserInfo
    # Performs authorization on user/collection

    def delete_collection(self, collection_id: str, user_info: UserInfo) -> None:
        pass

    # Update_collection
    # Replaces the current update_collection
    # Accepts:
    # Collection_id
    # UserInfo
    # A dataclass with the body to be updated
    # Should validate the body
    # Should handle DOI updates (re-use the existing logic with minimal refactors)
    # Can either return nothing or the metadata of the updated collection

    # TODO: body should be a dataclass?
    def update_collection_version(self, version_id: str, body: dict, user_info: UserInfo) -> None:
        pass

    # Create_collection_version
    # Replaces the current post_collection_revision
    # Accepts:
    # Collection_id
    # UserInfo
    # Performs authorization on the collection
    # Since revision logic is database specific, it delegates to the underlying layer
    # Returns a handle to the revised collection (either id or the full collection metadata)

    def create_collection_version(self, collection_id: str, user_info: UserInfo) -> CollectionVersion:
        pass

    def delete_collection_version(self, version_id: str, user_info: UserInfo) -> None:
        pass


    # Publish_collection
    # Replaces post (in publish.py)
    # Accepts:
    # Collection_id
    # UserInfo
    # Performs validation to make sure that the collection can be published
    # [Currently] accepts data_submission_policy_version: what is this for?
    # [Currently] triggers Cloudfront invalidation for the index endpoints. This should arguably NOT be done here but by the API layer
    # Since revision logic is database specific, it delegates to the underlying layer

    def publish_collection_version(self, version_id: str, user_info: UserInfo) -> None:
        pass

    # Ingest_dataset
    # Replaces the existing Upload_from_link
    # Potentially, also replaces relink (I am not sure why they are 2 separate functions)
    # Accepts:
    # Collection_id
    # UserInfo
    # URL of the uploadable dataset
    # [Optional] a dataset_id to be replaced
    # This is one of the most complex functions. Other than the database provider, It will need two additional providers:
    # StepFunctionProvider (to call the SFN that triggers the upload)
    # DropboxProvider to interface with Dropbox (could be more generic: RemoteFileProvider?)
    # Should handle exceptions from all providers:
    # Should only raise custom exceptions

    def ingest_dataset(self, collection_version_id: str, url: str, existing_dataset_version_id: Optional[str], user_info: UserInfo) -> str:
        pass

    # Get_all_datasets
    # Replaces get_dataset_index

    def get_all_datasets(self) -> Iterable[DatasetVersion]:
        pass

    # Delete_dataset
    # Replaces delete_dataset

    def delete_dataset(self, dataset_version_id: str) -> None:
        pass


    # get_dataset_assets
    # Replaces get_dataset_assets

    def get_dataset_artifacts(self, dataset_id: str) -> Iterable[DatasetArtifact]:
        pass


    # Download_dataset_asset
    # Replaces post_dataset_asset


    def download_dataset_asset(self, dataset_id: str) -> str:
        pass

    # Get_dataset_status
    # Replaces get_status

    def get_dataset_status(self, dataset_id: str) -> DatasetStatus:
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