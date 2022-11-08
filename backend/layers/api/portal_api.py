from typing import Optional
from unittest.mock import Mock

from flask import jsonify, make_response
from backend.corpora.common.utils.authorization_checks import is_user_owner_or_allowed
from backend.corpora.common.utils.http_exceptions import ForbiddenHTTPException
from backend.layers.auth.user_info import UserInfo
from backend.layers.business.business import BusinessLogic

from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.business.entities import CollectionQueryFilter
from backend.layers.common.entities import CollectionId, CollectionVersionId, DatasetArtifact, DatasetStatus, DatasetVersion, Link

from backend.corpora.common.utils import authorization_checks as auth
import itertools

from backend.layers.thirdparty.crossref_provider import CrossrefProviderInterface
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.step_function_provider import StepFunctionProviderInterface
from backend.layers.thirdparty.uri_provider import UriProviderInterface
from tests.unit.backend.layers.persistence.persistence_mock import DatabaseProviderMock


class PortalApi:

    business_logic: BusinessLogicInterface

    def __init__(self, business_logic: BusinessLogic) -> None:
        self.business_logic = business_logic

    def get_collections_list(self, from_date: int = None, to_date: int = None, token_info: Optional[dict] = None):
        """
        Returns all collections that are either published or belong to the user.
        `from_date` and `to_date` are deprecated parameters and should not be used.
        If there is no token_info, only published collections should be returned
        """

        all_published_collections = self.business_logic.get_collections(CollectionQueryFilter(is_published=True))

        user_info = UserInfo(token_info) # TODO: ideally, connexion should already return a UserInfo object

        if user_info.is_none():
            all_owned_collections = []
        elif user_info.is_super_curator():
            all_owned_collections = self.business_logic.get_collections(CollectionQueryFilter(is_published=False))
        else:
            all_owned_collections = self.business_logic.get_collections(CollectionQueryFilter(is_published=False, owner=user_info.user_id()))

        collections = []
        for c in itertools.chain(all_published_collections, all_owned_collections):
            collections.append({
                "id": c.collection_id.id,
                "visibility": "PRIVATE" if c.published_at is None else "PUBLIC",
                "owner": c.owner, # TODO: looks like this isn't returned right now
                "created_at": 12345, # TODO
                "revision_of": "TODO", # TODO: looks like this isn't returned right now
            })            

        result = {"collections": collections}
        return make_response(jsonify(result), 200)

    def _dataset_processing_status_to_response(self, status: DatasetStatus, dataset_id: str):
        return {
            "created_at": 1234,
            "cxg_status": status.cxg_status,
            "dataset_id": dataset_id,
            "h5ad_status": status.h5ad_status,
            "id": "TODO", # TODO can we purge?
            "processing_status": status.processing_status,
            "rds_status": status.rds_status,
            "updated_at": 1234,
            "upload_progress": 1234,
            "upload_status": status.upload_status,
            "validation_status": status.validation_status,
        }

    def _link_to_response(self, link: Link):
        return {
            "link_name": link.name,
            "link_type": link.type,
            "link_url": link.uri,
        }

    def _dataset_asset_to_response(self, dataset_artifact: DatasetArtifact, dataset_id: str):
        return {
            "created_at": 1234,
            "dataset_id": dataset_id,
            "filename": "TODO",
            "filetype": dataset_artifact.type,
            "id": dataset_artifact.id,
            "s3_uri": dataset_artifact.uri,
            "updated_at": 1234,
            "user_submitted": True,
        }

    def _dataset_to_response(self, dataset: DatasetVersion):
        return {
            "assay": dataset.metadata.assay,
            "batch_condition": dataset.metadata.batch_condition,
            "cell_count": dataset.metadata.cell_count,
            "cell_type": dataset.metadata.cell_type,
            "collection_id": "TODO", # TODO
            "created_at": 1234, # TODO
            "dataset_assets": [self._dataset_asset_to_response(a, dataset.dataset_id) for a in dataset.artifacts],
            "dataset_deployments": dataset.metadata.explorer_url,
            "development_stage": dataset.metadata.deployment_stage,
            "disease": dataset.metadata.disease,
            "donor_id": dataset.metadata.donor_id,
            "id": dataset.id,
            "is_primary_data": dataset.metadata.is_primary_data,
            "is_valid": True, # why do we have this
            "mean_genes_per_cell": dataset.metadata.mean_genes_per_cell,
            "name": dataset.metadata.name,
            "organism": dataset.metadata.organism,
            "processing_status": self._dataset_processing_status_to_response(dataset.status, dataset.dataset_id),
            "published": True,
            "published_at": 1234,
            "revision": 1234,
            "schema_version": "3.0.0",
            "self_reported_ethnicity": dataset.metadata.self_reported_ethnicity,
            "sex": dataset.metadata.sex,
            "suspension_type": dataset.metadata.suspension_type,
            "tissue": dataset.metadata.tissue,
            "tombstone": False,
            "updated_at": 1234,
            "x_approximate_distribution": dataset.metadata.x_approximate_distribution,
        }

    def get_collection_details(self, collection_id: str, token_info: dict):
        """
        Retrieves the collection information. Will look up for a published collection first,
        and then looks up for a collection version
        """
        # TODO: this logic might belong to the business layer?
        version = self.business_logic.get_published_collection_version(CollectionId(collection_id))
        if version is None:
            version = self.business_logic.get_collection_version(CollectionVersionId(collection_id))

        print(version)
        
        if version is None:
            raise ForbiddenHTTPException() # TODO: maybe remake this collection

        # TODO: handle tombstoning

        # TODO: handle access type as follows
        #result["access_type"] = "WRITE" if is_user_owner_or_allowed(token_info, collection.owner) else "READ"

        response = {
            "access_type": "TODO", # TODO
            "contact_email": "example@gmail.com",
            "contact_name": "Smoke Test User",
            "created_at": 1234,
            "curator_name": "", # TODO
            "data_submission_policy_version": "1.0", # TODO
            "datasets": [self._dataset_to_response(d) for d in version.datasets],
            "description": version.metadata.description,
            "id": version.collection_id,
            "links": [self._link_to_response(l) for l in version.metadata.links],
            "name": version.metadata.name,
            "published_at": 1234,
            "updated_at": 1234,
            "visibility": "PUBLIC" if version.published_at is not None else "PRIVATE",
        }
 
        return make_response(jsonify(response), 200)

        
        # db_session = g.db_session
        # collection = get_collection_else_forbidden(db_session, collection_id, include_tombstones=True)
        # if collection.tombstone:
        #     result = ""
        #     response = 410
        # else:
        #     get_tombstone_datasets = (
        #         is_user_owner_or_allowed(token_info, collection.owner)
        #         and collection.visibility == CollectionVisibility.PRIVATE
        #     )
        #     result = collection.reshape_for_api(get_tombstone_datasets)
        #     response = 200
        #     result["access_type"] = "WRITE" if is_user_owner_or_allowed(token_info, collection.owner) else "READ"
        # return make_response(jsonify(result), response)

    def get_collections_index(self):
        pass

    def post_collection_revision(self, collection_id: str, token_info: dict):
        pass

    def create_collection(self, body: dict, user: str):
        pass

    def delete_collection(self, collection_id: str, token_info: dict):
        pass

    def update_collection(self, collection_id: str, body: dict, token_info: dict):
        pass


    def post(self, collection_id: str, body: object, token_info: dict): # publish
        pass

    def link(self, collection_id: str, body: dict, token_info: dict):
        pass

    def relink(self, collection_id: str, body: dict, token_info: dict):
        pass

    def upload_from_link(self, collection_id: str, token_info: dict, url: str, dataset_id: str = None):
        pass



    def post_dataset_asset(self, dataset_id: str, asset_id: str):
        pass

    def get_dataset_assets(self, dataset_id: str):
        pass

    def get_status(self, dataset_id: str, token_info: dict):
        pass

    def get_datasets_index(self):
        pass

    def delete_dataset(self, dataset_id: str, token_info: dict):
        pass

    def get_dataset_identifiers(self, url: str):
        pass

    def post_dataset_gene_sets(self, dataset_id: str, body: object, token_info: dict):
        pass
