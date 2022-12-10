from dataclasses import asdict

from flask import Response, jsonify, make_response

from backend.common.corpora_config import CorporaConfig
from backend.common.corpora_orm import ProjectLinkType
from backend.common.utils.http_exceptions import (
    InvalidParametersHTTPException,
    MethodNotAllowedException,
)
from backend.layers.api.router import get_business_logic
from backend.layers.auth.user_info import UserInfo
from backend.layers.business.entities import CollectionMetadataUpdate
from backend.layers.business.exceptions import CollectionUpdateException
from backend.layers.common import doi
from backend.layers.common.entities import Link
from backend.portal.api.curation.v1.curation.collections.common import (
    allowed_dataset_asset_types,
    extract_doi_from_links,
    get_collection_level_processing_status,
    get_infered_collection_version_else_forbidden,
    get_visibility,
    is_owner_or_allowed_else_forbidden,
    is_primary_data_mapping,
)


def delete(collection_id: str, token_info: dict) -> Response:
    user_info = UserInfo(token_info)
    collection_version = get_infered_collection_version_else_forbidden(collection_id)
    is_owner_or_allowed_else_forbidden(collection_version, user_info)
    if collection_version.published_at:
        raise MethodNotAllowedException(detail="Cannot delete a published collection through API.")
    else:
        get_business_logic().delete_collection_version(collection_version.version_id)
    return make_response("", 204)


def get(collection_id: str, token_info: dict) -> Response:
    collection_version = get_infered_collection_version_else_forbidden(collection_id)
    user_info = UserInfo(token_info)

    # get collectoin attributes based on published status
    if collection_version.published_at is None:
        # Unpublished
        resp_collection_id = collection_version.version_id
        revision_of = collection_version.collection_id.id
        revising_in = None
    else:
        # Published
        resp_collection_id = collection_version.collection_id
        revision_of = None
        if not user_info.is_user_owner_or_allowed(collection_version.owner):
            _revising_in = None
        else:
            _revising_in = get_business_logic().get_unplublished_collection_version_from_canonical(resp_collection_id)
        revising_in = _revising_in.version_id.id if _revising_in else None

    # get collection dataset attributes
    response_datasets = []
    for dataset in collection_version.datasets:
        ds = asdict(dataset.metadata)

        # get dataset asset attributes
        assets = []
        for artifact in dataset.artifacts:
            if artifact.type in allowed_dataset_asset_types:
                assets.append(dict(filetype=artifact.type, filename=artifact.uri.split("/")[-1]))
        ds["dataset_assets"] = assets
        ds["processing_status_detail"] = dataset.status.validation_message
        ds["revised_at"] = dataset.canonical_dataset.revised_at
        ds["is_primary_data"] = is_primary_data_mapping.get(ds.pop("is_primary_data"), [])
        response_datasets.append(ds)

    # build response
    processing_status = get_collection_level_processing_status(collection_version.datasets)
    revised_at = get_business_logic().get_published_collection_version(collection_version.collection_id).published_at
    doi, links = extract_doi_from_links(collection_version.metadata.links)
    response = dict(
        collection_url=f"{CorporaConfig().collections_base_url}/collections/{collection_id}",
        contact_email=collection_version.metadata.contact_email,
        contact_name=collection_version.metadata.contact_name,
        created_at=collection_version.created_at,
        curator_name=collection_version.owner,
        datasets=response_datasets,
        description=collection_version.metadata.description,
        dio=doi,
        id=collection_id,
        links=links,
        name=collection_version.metadata.name,
        processing_status=processing_status,
        published_at=collection_version.canonical_collection.originally_published_at,
        publisher_metadata=collection_version.publisher_metadata,
        revised_at=revised_at,
        revising_in=revising_in,
        revision_of=revision_of,
        visibility=get_visibility(collection_version),
    )
    return jsonify(response)


def patch(collection_id: str, body: dict, token_info: dict) -> Response:
    user_info = UserInfo(token_info)
    collection_version = get_infered_collection_version_else_forbidden(collection_id)
    is_owner_or_allowed_else_forbidden(collection_version, user_info)
    if collection_version.published_at:
        raise MethodNotAllowedException(
            detail="Directly editing a public Collection is not allowed; you must create a revision."
        )

    # Build CollectionMetadataUpdate object
    body["links"] = [Link(link.get("link_name"), link["link_type"], link["link_url"]) for link in body.get("links", [])]
    collection_metadata = CollectionMetadataUpdate(**body)

    # Update the collection
    errors = []
    if doi_url := body.get("doi"):
        if doi_url := doi.curation_get_normalized_doi_url(doi_url, errors):
            links = body.get("links", [])
            links.append({"link_type": ProjectLinkType.DOI.name, "link_url": doi_url})
            body["links"] = links
    try:
        get_business_logic().update_collection_version(collection_version.version_id, collection_metadata)
    except CollectionUpdateException as ex:
        errors.extend(ex.errors)
    if errors:
        raise InvalidParametersHTTPException(ext=dict(invalid_parameters=errors))

    # Make Response
    updated_collection_version = get_business_logic().get_collection_version(collection_version.version_id)

    metadata = asdict(updated_collection_version.metadata)
    metadata.pop("links")

    doi_url, links = extract_doi_from_links(updated_collection_version.metadata.links)

    response = dict(**metadata, publisher_metadata=updated_collection_version.publisher_metadata, links=links, doi=doi)
    return jsonify(response)
