import itertools
from typing import List, Optional, Tuple
from urllib.parse import urlparse

from flask import Response, jsonify, make_response

from backend.common.utils.http_exceptions import (
    ConflictException,
    ForbiddenHTTPException,
    InvalidParametersHTTPException,
    MethodNotAllowedException,
    NotFoundHTTPException,
    ServerErrorHTTPException,
    TooLargeHTTPException,
)
from backend.common.utils.ontology_mappings.ontology_map_loader import ontology_mappings
from backend.layers.api.common import ApiCommon
from backend.layers.api.enrichment import enrich_dataset_with_ancestors
from backend.layers.auth.user_info import UserInfo
from backend.layers.business.business import BusinessLogic
from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.business.entities import CollectionMetadataUpdate, CollectionQueryFilter
from backend.layers.business.exceptions import (
    ArtifactNotFoundException,
    CollectionCreationException,
    CollectionIsPublishedException,
    CollectionNotFoundException,
    CollectionPublishException,
    CollectionUpdateException,
    CollectionVersionException,
    DatasetInWrongStatusException,
    DatasetNotFoundException,
    InvalidURIException,
    MaxFileSizeExceededException,
)
from backend.layers.common import doi
from backend.layers.common.entities import (
    CollectionId,
    CollectionMetadata,
    CollectionVersion,
    CollectionVersionId,
    DatasetArtifact,
    DatasetArtifactId,
    DatasetStatus,
    DatasetVersion,
    DatasetVersionId,
    Link,
    OntologyTermId,
)


class PortalApi(ApiCommon):
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

        user_info = UserInfo(token_info)  # TODO: ideally, connexion should already return a UserInfo object

        if user_info.is_none():
            all_owned_collections = []
        elif user_info.is_super_curator():
            all_owned_collections = self.business_logic.get_collections(CollectionQueryFilter(is_published=False))
        else:
            all_owned_collections = self.business_logic.get_collections(
                CollectionQueryFilter(is_published=False, owner=user_info.user_id)
            )

        collections = []
        for c in itertools.chain(all_published_collections, all_owned_collections):
            collection = {
                "id": c.version_id.id if c.published_at is None else c.collection_id.id,
                "visibility": "PRIVATE" if c.published_at is None else "PUBLIC",
                "owner": c.owner,
                "created_at": c.created_at,
            }
            if c.published_at is None:
                collection["revision_of"] = c.collection_id.id
            collections.append(collection)

        result = {"collections": collections}
        return make_response(jsonify(result), 200)

    def _dataset_processing_status_to_response(self, status: DatasetStatus, dataset_id: str):
        return {
            "created_at": 0,  # NA
            "cxg_status": status.cxg_status or "NA",
            "dataset_id": dataset_id,
            "h5ad_status": status.h5ad_status or "NA",
            "id": "NA",  # TODO can we purge?
            "processing_status": status.processing_status or "NA",
            "rds_status": status.rds_status or "NA",
            "updated_at": 0,  # NA
            "upload_progress": 1,  # No longer supported - always return 1
            "upload_status": status.upload_status or "NA",
            "validation_status": status.validation_status or "NA",
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
            "created_at": 0,
            "dataset_id": dataset_id,
            "filename": "TODO",  # TODO: might need to get it from the url
            "filetype": dataset_artifact.type.upper(),
            "id": dataset_artifact.id.id,
            "s3_uri": dataset_artifact.uri,
            "updated_at": 0,
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

    # Note: `metadata` can be none while the dataset is uploading
    def _dataset_to_response(self, dataset: DatasetVersion):
        return self.remove_none(
            {
                "assay": None
                if dataset.metadata is None
                else self._ontology_term_ids_to_response(dataset.metadata.assay),
                "batch_condition": None if dataset.metadata is None else dataset.metadata.batch_condition,
                "cell_count": None if dataset.metadata is None else dataset.metadata.cell_count,
                "cell_type": None
                if dataset.metadata is None
                else self._ontology_term_ids_to_response(dataset.metadata.cell_type),
                "collection_id": dataset.collection_id.id,
                "created_at": dataset.created_at,
                "dataset_assets": [
                    self._dataset_asset_to_response(a, dataset.version_id.id) for a in dataset.artifacts
                ],
                "dataset_deployments": [{"url": "TODO"}],  # TODO: dataset.metadata.explorer_url,
                "development_stage": None
                if dataset.metadata is None
                else self._ontology_term_ids_to_response(dataset.metadata.development_stage),
                "disease": None
                if dataset.metadata is None
                else self._ontology_term_ids_to_response(dataset.metadata.disease),
                "donor_id": None if dataset.metadata is None else dataset.metadata.donor_id,
                "id": dataset.version_id.id,
                "is_primary_data": None if dataset.metadata is None else dataset.metadata.is_primary_data,
                "is_valid": True,  # why do we have this
                "mean_genes_per_cell": None if dataset.metadata is None else dataset.metadata.mean_genes_per_cell,
                "name": "" if dataset.metadata is None else dataset.metadata.name,
                "organism": None
                if dataset.metadata is None
                else self._ontology_term_ids_to_response(dataset.metadata.organism),
                "processing_status": self._dataset_processing_status_to_response(dataset.status, dataset.version_id.id),
                "published": True,  # TODO
                "published_at": dataset.canonical_dataset.published_at,
                "revision": 0,  # TODO this is the progressive revision number. I don't think we'll need this
                "schema_version": None if dataset.metadata is None else dataset.metadata.schema_version,
                "self_reported_ethnicity": None
                if dataset.metadata is None
                else self._ontology_term_ids_to_response(dataset.metadata.self_reported_ethnicity),
                "sex": None if dataset.metadata is None else self._ontology_term_ids_to_response(dataset.metadata.sex),
                "suspension_type": None if dataset.metadata is None else dataset.metadata.suspension_type,
                "tissue": None
                if dataset.metadata is None
                else self._ontology_term_ids_to_response(dataset.metadata.tissue),
                "tombstone": False,  # TODO
                "updated_at": dataset.created_at,  # Legacy: datasets can't be updated anymore
                "x_approximate_distribution": None
                if dataset.metadata is None
                else dataset.metadata.x_approximate_distribution,
            }
        )

    def _collection_to_response(self, collection: CollectionVersion, access_type: str):
        collection_id = collection.collection_id.id if collection.published_at is not None else collection.version_id.id
        return self.remove_none(
            {
                "access_type": access_type,
                "contact_email": collection.metadata.contact_email,
                "contact_name": collection.metadata.contact_name,
                "created_at": collection.created_at,
                "curator_name": "",  # TODO
                "data_submission_policy_version": "1.0",  # TODO
                "datasets": [self._dataset_to_response(ds) for ds in collection.datasets],
                "description": collection.metadata.description,
                "id": collection_id,
                "links": [self._link_to_response(link) for link in collection.metadata.links],
                "name": collection.metadata.name,
                "published_at": collection.published_at,
                "publisher_metadata": collection.publisher_metadata,  # TODO: convert
                "updated_at": collection.published_at or collection.created_at,
                "visibility": "PUBLIC" if collection.published_at is not None else "PRIVATE",
            }
        )

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
            raise ForbiddenHTTPException()  # TODO: maybe remake this exception

        # TODO: handle tombstoning

        user_info = UserInfo(token_info)
        print(token_info, user_info.user_id, version.owner)
        access_type = "WRITE" if user_info.is_user_owner_or_allowed(version.owner) else "READ"

        response = self._collection_to_response(version, access_type)

        return make_response(jsonify(response), 200)

    def post_collection_revision(self, collection_id: str, token_info: dict):
        """
        Creates a collection revision
        """

        published_collection = self.business_logic.get_published_collection_version(CollectionId(collection_id))

        if published_collection is None:
            raise ForbiddenHTTPException(f"Collection {collection_id} does not exist")

        if not UserInfo(token_info).is_user_owner_or_allowed(published_collection.owner):
            raise ForbiddenHTTPException("Unauthorized")

        try:
            version = self.business_logic.create_collection_version(CollectionId(collection_id))
        except CollectionVersionException:
            raise ForbiddenHTTPException("Another revision is already in progress")

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
            raise InvalidParametersHTTPException(detail=errors)  # TODO: rewrite this exception?

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
                "published_at": collection.canonical_collection.originally_published_at,
                "revised_at": collection.published_at,
            }

            if collection.publisher_metadata is not None:
                transformed_collection["publisher_metadata"] = self._publisher_metadata_to_response(
                    collection.publisher_metadata
                )

            response.append(transformed_collection)

        return make_response(jsonify(response), 200)

    def delete_collection(self, collection_id: str, token_info: dict):
        """
        Deletes a collection version from the persistence store.
        """
        version = self.business_logic.get_collection_version(CollectionVersionId(collection_id))
        if not UserInfo(token_info).is_user_owner_or_allowed(version.owner):
            raise ForbiddenHTTPException()

        self.business_logic.delete_collection_version(CollectionVersionId(collection_id))

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

    def publish_post(self, collection_id: str, body: object, token_info: dict):
        """
        Publishes a collection
        """

        version = self.business_logic.get_collection_version(CollectionVersionId(collection_id))
        if version is None or not UserInfo(token_info).is_user_owner_or_allowed(version.owner):
            raise ForbiddenHTTPException()

        try:
            self.business_logic.publish_collection_version(CollectionVersionId(collection_id))
        except CollectionPublishException:
            raise ConflictException(detail="The collection must have a least one dataset.")

        return make_response({"collection_id": version.collection_id.id, "visibility": "PUBLIC"}, 202)

    def upload_from_link(self, collection_id: str, token_info: dict, url: str, dataset_id: str = None):

        version = self.business_logic.get_collection_version(CollectionVersionId(collection_id))
        if version is None or not UserInfo(token_info).is_user_owner_or_allowed(version.owner):
            raise ForbiddenHTTPException()

        try:
            dataset_version_id, _ = self.business_logic.ingest_dataset(
                CollectionVersionId(collection_id),
                url,
                None if dataset_id is None else DatasetVersionId(dataset_id),
            )
            return dataset_version_id
        except CollectionNotFoundException:
            raise ForbiddenHTTPException()
        except CollectionIsPublishedException:
            raise ForbiddenHTTPException()
        except DatasetNotFoundException:
            raise NotFoundHTTPException()
        except InvalidURIException:
            raise InvalidParametersHTTPException(detail="The dropbox shared link is invalid.")
        except MaxFileSizeExceededException:
            raise TooLargeHTTPException()
        except DatasetInWrongStatusException:
            raise MethodNotAllowedException(
                detail="Submission failed. A dataset cannot be updated while a previous update for the same dataset "
                "is in progress. Please cancel the current submission by deleting the dataset, or wait until "
                "the submission has finished processing."
            )

    # TODO: those two methods should probably be collapsed into one
    # TODO: not quite sure what's the difference between url and link - investigate
    def upload_link(self, collection_id: str, body: dict, token_info: dict):
        dataset_id = self.upload_from_link(collection_id, token_info, body["url"])
        return make_response({"dataset_id": dataset_id.id}, 202)

    def upload_relink(self, collection_id: str, body: dict, token_info: dict):
        dataset_id = self.upload_from_link(
            collection_id,
            token_info,
            body.get("url", body.get("link")),
            body.get("id"),
        )
        return make_response({"dataset_id": dataset_id.id}, 202)

    def post_dataset_asset(self, dataset_id: str, asset_id: str):
        """
        Requests to download a dataset asset, by generating a presigned_url.
        """

        version = self.business_logic.get_dataset_version(DatasetVersionId(dataset_id))
        if version is None:
            raise NotFoundHTTPException(detail=f"'dataset/{dataset_id}' not found.")

        try:
            download_data = self.business_logic.get_dataset_artifact_download_data(
                DatasetVersionId(dataset_id), DatasetArtifactId(asset_id)
            )
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
            artifacts.append({"id": artifact.id.id, "type": artifact.type, "s3_uri": artifact.uri})
        response = {"assets": artifacts}

        return make_response(jsonify(response), 200)

    def get_status(self, dataset_id: str, token_info: dict):
        """
        Returns the processing status of a dataset.
        Raises Unauthorized if the dataset does not exist or if the user is not authorized for it
        """
        version, _ = self._assert_dataset_has_right_owner(DatasetVersionId(dataset_id), UserInfo(token_info))

        response = {
            "cxg_status": version.status.cxg_status or "NA",
            "rds_status": version.status.rds_status or "NA",
            "h5ad_status": version.status.h5ad_status or "NA",
            "processing_status": version.status.processing_status or "NA",
            "dataset_id": dataset_id,
            "id": "NA",
            "upload_progress": 1,
            "upload_status": version.status.upload_status or "NA",
            "validation_status": version.status.validation_status or "NA",
        }

        return make_response(response, 200)

    def get_datasets_index(self):
        """
        Returns a list of all the datasets that currently belong to a published and active collection
        """

        response = []
        for dataset in self.business_logic.get_all_published_datasets():
            payload = self._dataset_to_response(dataset)
            enrich_dataset_with_ancestors(
                payload, "development_stage", ontology_mappings.development_stage_ontology_mapping
            )
            enrich_dataset_with_ancestors(payload, "tissue", ontology_mappings.tissue_ontology_mapping)
            enrich_dataset_with_ancestors(payload, "cell_type", ontology_mappings.cell_type_ontology_mapping)
            response.append(payload)

        return make_response(jsonify(response), 200)

    def delete_dataset(self, dataset_id: str, token_info: dict):
        """
        Deletes a dataset version.
        """
        dataset_version, collection_version = self._assert_dataset_has_right_owner(
            DatasetVersionId(dataset_id), UserInfo(token_info)
        )
        if dataset_version.version_id not in [v.version_id for v in collection_version.datasets]:
            raise ForbiddenHTTPException(f"Dataset {dataset_id} does not belong to a collection")

        try:
            self.business_logic.remove_dataset_version(collection_version.version_id, DatasetVersionId(dataset_id))
        except CollectionUpdateException:
            raise MethodNotAllowedException(detail="Cannot delete a public Dataset")
        return Response(status=202)

    def get_dataset_identifiers(self, url: str):
        """
        a.k.a. the meta endpoint
        """
        try:
            path = urlparse(url).path
            id = [segment for segment in path.split("/") if segment][-1].strip(".cxg")
        except Exception:
            raise ServerErrorHTTPException("Cannot parse URL")

        # TODO: if we require it, add double id lookup
        dataset = self.business_logic.get_dataset_version(DatasetVersionId(id))
        if dataset is None:
            raise NotFoundHTTPException()

        collection = self.business_logic.get_collection_version_from_canonical(dataset.collection_id)
        if collection is None:  # orphaned datasets
            raise NotFoundHTTPException()

        # Retrieves the URI of the cxg artifact
        s3_uri = next(a.uri for a in dataset.artifacts if a.type == "CXG")

        dataset_identifiers = {
            "s3_uri": s3_uri,
            "dataset_id": dataset.dataset_id.id,
            "collection_id": dataset.collection_id.id,
            "collection_visibility": "PUBLIC" if collection.published_at is not None else "PRIVATE",
            "tombstoned": False,  # No longer applicable
        }
        return make_response(jsonify(dataset_identifiers), 200)
