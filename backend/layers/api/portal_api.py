from typing import List, Optional
from unittest.mock import Mock

from flask import jsonify, make_response
from backend.corpora.common.utils.authorization_checks import is_user_owner_or_allowed
from backend.corpora.common.utils.http_exceptions import ForbiddenHTTPException, InvalidParametersHTTPException, NotFoundHTTPException, ServerErrorHTTPException
from backend.layers.api.enrichment import enrich_dataset_with_ancestors
from backend.layers.auth.user_info import UserInfo
from backend.layers.business.business import BusinessLogic

from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.business.entities import CollectionMetadataUpdate, CollectionQueryFilter
from backend.layers.business.exceptions import ArtifactNotFoundException, CollectionCreationException
from backend.layers.common.entities import CollectionId, CollectionMetadata, CollectionVersion, CollectionVersionId, DatasetArtifact, DatasetId, DatasetStatus, DatasetVersion, DatasetVersionId, Link, OntologyTermId

from backend.corpora.common.utils import authorization_checks as auth
from backend.corpora.common.utils.ontology_mappings.ontology_map_loader import ontology_mappings
import itertools

from backend.layers.common import doi

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
                # "id": c.collection_id.id,
                "id": c.version_id.id if c.published_at is None else c.collection_id.id,
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

    # TODO: use remove_none
    def _link_to_response(self, link: Link):
        response = {
            "link_type": link.type,
            "link_url": link.uri,
        }
        if link.name is not None:
            response["link_name"] = link.name
        return response


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

    def _ontology_term_id_to_response(self, ontology_term_id: OntologyTermId):
        return {
            "label": ontology_term_id.label,
            "ontology_term_id": ontology_term_id.ontology_term_id,
        }

    def _ontology_term_ids_to_response(self, ontology_term_ids: List[OntologyTermId]):
        return [self._ontology_term_id_to_response(otid) for otid in ontology_term_ids]

    def remove_none(self, body: dict):
        return {k: v for k, v in body.items() if v is not None}

    def _dataset_to_response(self, dataset: DatasetVersion):
        return {
            "assay": self._ontology_term_ids_to_response(dataset.metadata.assay),
            "batch_condition": dataset.metadata.batch_condition,
            "cell_count": dataset.metadata.cell_count,
            "cell_type": self._ontology_term_ids_to_response(dataset.metadata.cell_type),
            "collection_id": "TODO", # TODO
            "created_at": 1234, # TODO
            "dataset_assets": [self._dataset_asset_to_response(a, dataset.dataset_id.id) for a in dataset.artifacts],
            "dataset_deployments": [{"url": "TODO"}], # TODO: dataset.metadata.explorer_url,
            "development_stage": self._ontology_term_ids_to_response(dataset.metadata.development_stage),
            "disease": self._ontology_term_ids_to_response(dataset.metadata.disease),
            "donor_id": dataset.metadata.donor_id,
            "id": dataset.dataset_id.id,
            "is_primary_data": dataset.metadata.is_primary_data,
            "is_valid": True, # why do we have this
            "mean_genes_per_cell": dataset.metadata.mean_genes_per_cell,
            "name": dataset.metadata.name,
            "organism": self._ontology_term_ids_to_response(dataset.metadata.organism),
            "processing_status": self._dataset_processing_status_to_response(dataset.status, dataset.dataset_id.id),
            "published": True,
            "published_at": 1234,
            "revision": 1234,
            "schema_version": "3.0.0",
            "self_reported_ethnicity": self._ontology_term_ids_to_response(dataset.metadata.self_reported_ethnicity),
            "sex": self._ontology_term_ids_to_response(dataset.metadata.sex),
            "suspension_type": dataset.metadata.suspension_type,
            "tissue": self._ontology_term_ids_to_response(dataset.metadata.tissue),
            "tombstone": False,
            "updated_at": 1234,
            "x_approximate_distribution": dataset.metadata.x_approximate_distribution,
        }

    def _collection_to_response(self, collection: CollectionVersion, access_type: str):
        collection_id = collection.collection_id.id if collection.published_at is not None else collection.version_id.id
        return self.remove_none({
            "access_type": access_type,
            "contact_email": collection.metadata.contact_email,
            "contact_name": collection.metadata.contact_name,
            "created_at": 1234,
            "curator_name": "", # TODO
            "data_submission_policy_version": "1.0", # TODO
            "datasets": [self._dataset_to_response(d) for d in collection.datasets],
            "description": collection.metadata.description,
            "id": collection_id,
            "links": [self._link_to_response(l) for l in collection.metadata.links],
            "name": collection.metadata.name,
            "published_at": 1234,
            "publisher_metadata": collection.publisher_metadata, # TODO: convert
            "updated_at": 1234,
            "visibility": "PUBLIC" if collection.published_at is not None else "PRIVATE",
        })

    def get_collection_details(self, collection_id: str, token_info: dict):
        """
        Retrieves the collection information. Will look up for a published collection first,
        and then looks up for a collection version
        """
        # TODO: this logic might belong to the business layer?
        version = self.business_logic.get_published_collection_version(CollectionId(collection_id))
        if version is None:
            version = self.business_logic.get_collection_version(CollectionVersionId(collection_id))
        if version is None:
            raise ForbiddenHTTPException() # TODO: maybe remake this exception

        # TODO: handle tombstoning

        user_info = UserInfo(token_info)
        print(token_info, user_info.user_id(), version.owner)
        access_type = "WRITE" if user_info.is_user_owner_or_allowed(version.owner) else "READ"

        response = self._collection_to_response(version, access_type)
 
        return make_response(jsonify(response), 200)


    def post_collection_revision(self, collection_id: str, token_info: dict):
        version = self.business_logic.create_collection_version(CollectionId(collection_id))
        response = self._collection_to_response(version, "WRITE")

        return make_response(response, 201)

    def _link_from_request(self, body: dict):
        return Link(
            body.get("link_name"),
            body["link_type"],
            body["link_url"],
        )

    # TODO: why do we have `user` and not `token_info`? This seems weird
    def create_collection(self, body: dict, user: str):
        """
        Creates a collection. Will also perform DOI normalization: if the DOI is specified in `links`
        as a CURIE (i.e., without the https://doi.org prefix), it will be normalized.
        All exceptions are caught and raised as an InvalidParametersHTTPException.
        """

        errors = []
        doi_url = None
        if doi_node := doi.get_doi_link_node(body, errors):
            if doi_url := doi.portal_get_normalized_doi_url(doi_node, errors):
                doi_node["link_url"] = doi_url

        if errors:
            raise InvalidParametersHTTPException(detail=errors) # TODO: rewrite this exception?

        metadata = CollectionMetadata(
            body["name"],
            body["description"],
            body["contact_name"],
            body["contact_email"],
            [self._link_from_request(node) for node in body.get("links", [])],
        )

        try:
            version = self.business_logic.create_collection(user, metadata)
        except CollectionCreationException as ex:
            raise InvalidParametersHTTPException(detail=ex.errors)

        return make_response(jsonify({"collection_id": version.version_id.id}), 201)

    # TODO: we should use a dataclass here
    def _publisher_metadata_to_response(self, publisher_metadata: dict):
        return publisher_metadata

    def get_collection_index(self):
        """
        Returns a list of collections that are published and active. 
        Also returns a subset of fields and not datasets.
        """
        collections = self.business_logic.get_collections(CollectionQueryFilter(is_published=True))
        response = []

        for collection in collections:
            transformed_collection = {
                "id": collection.collection_id.id,
                "name": collection.metadata.name,
                "published_at": collection.published_at,
                "revised_at": 0, # TODO
            }

            if collection.publisher_metadata is not None:
                transformed_collection["publisher_metadata"] = self._publisher_metadata_to_response(collection.publisher_metadata)

            response.append(transformed_collection)

        return make_response(jsonify(response), 200)


    def delete_collection(self, collection_id: str, token_info: dict):
        pass

    def update_collection(self, collection_id: str, body: dict, token_info: dict):
        """
        Updates a collection
        """

        # Ensure that the version exists and the user is authorized to update it
        # TODO: this should be extracted to a method, I think
        version = self.business_logic.get_collection_version(CollectionVersionId(collection_id))
        if version is None or not UserInfo(token_info).is_user_owner_or_allowed(version.owner):
            raise ForbiddenHTTPException()

        payload = CollectionMetadataUpdate(
            body.get("name"),
            body.get("description"),
            body.get("contact_name"),
            body.get("contact_email"),
            [self._link_from_request(node) for node in body.get("links", [])],
        )

        self.business_logic.update_collection_version(CollectionVersionId(collection_id), payload)

        # Requires strong consistency w.r.t. the operation above - if not available, the update needs 
        # to be done in memory
        version = self.business_logic.get_collection_version(CollectionVersionId(collection_id))

        response = self._collection_to_response(version, "WRITE")
        return make_response(jsonify(response), 200)


    def post(self, collection_id: str, body: object, token_info: dict): # publish
        pass

    def link(self, collection_id: str, body: dict, token_info: dict):
        pass

    def relink(self, collection_id: str, body: dict, token_info: dict):
        pass

    def upload_from_link(self, collection_id: str, token_info: dict, url: str, dataset_id: str = None):
        pass



    def post_dataset_asset(self, dataset_id: str, asset_id: str):
        """
        Requests to download a dataset asset, by generating a presigned_url.
        This method will accept a canonical dataset_id
        """

        # TODO: if dataset not found, raise NotFoundHTTPException(detail=f"'dataset/{dataset_id}' not found.")
        try:
            download_data = self.business_logic.get_dataset_artifact_download_data(DatasetId(dataset_id), asset_id)
        except ArtifactNotFoundException:
            raise NotFoundHTTPException(detail=f"'dataset/{dataset_id}/asset/{asset_id}' not found.")

        if download_data.file_size is None:
            raise ServerErrorHTTPException()

        if download_data.presigned_url is None:
            raise ServerErrorHTTPException()

        response = {
            "dataset_id": dataset_id,
            "file_name": download_data.file_name,
            "file_size": download_data.file_size,
            "presigned_url": download_data.presigned_url,
        }

        return make_response(response, 200)
        

    def get_dataset_assets(self, dataset_id: str):
        """
        Returns a list of all the artifacts registered to a dataset.
        TODO: not sure where this is used and what the response should be
        """

        artifacts = []
        for artifact in self.business_logic.get_dataset_artifacts(DatasetVersionId(dataset_id)):
            artifacts.append({
                "id": artifact.id.id,
                "type": artifact.type,
                "s3_uri": artifact.uri
            })
        response = {"assets": artifacts}

        return make_response(jsonify(response), 200)

    def get_status(self, dataset_id: str, token_info: dict):

        # version = self.business_logic.get_dataset_version(DatasetVersionId(dataset_id))
        # collection_version = self.business_logic.get_collection_version(version)

        # if version is None:
            # return 

        # TODO: needs a review
        # TODO: this needs to access the collection that owns this dataset - but we don't have that key right now...

        response = {
            "cxg_status": status.cxg_status or "NA",
            "rds_status": status.rds_status or "NA",
            "h5ad_status": status.h5ad_status or "NA",
            "processing_status": status.processing_status or "NA",
            "dataset_id": dataset_id,
            "id": "NA",
            "upload_progress": 12345,
            "upload_status": status.upload_status or "NA",
            "validation_status": status.validation_status or "NA",
        }

        return make_response(response, 200)

    def get_datasets_index(self):

        response = []
        for dataset in self.business_logic.get_all_published_datasets():
            payload = self._dataset_to_response(dataset)
            enrich_dataset_with_ancestors(payload, "development_stage", ontology_mappings.development_stage_ontology_mapping)
            enrich_dataset_with_ancestors(payload, "tissue", ontology_mappings.tissue_ontology_mapping)
            enrich_dataset_with_ancestors(payload, "cell_type", ontology_mappings.cell_type_ontology_mapping)
            response.append(payload)

        return make_response(jsonify(response), 200)

    def delete_dataset(self, dataset_id: str, token_info: dict):
        pass

    def get_dataset_identifiers(self, url: str):
        pass

    def post_dataset_gene_sets(self, dataset_id: str, body: object, token_info: dict):
        pass
