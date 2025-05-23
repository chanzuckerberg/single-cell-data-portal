from dataclasses import asdict

from flask import Response, jsonify, make_response

import backend.common.doi as doi
from backend.common.utils.http_exceptions import (
    ForbiddenHTTPException,
    InvalidParametersHTTPException,
    MethodNotAllowedException,
)
from backend.curation.api.v1.curation.collections.common import (
    extract_doi_from_links,
    get_inferred_collection_version,
    is_owner_or_allowed_else_forbidden,
    reshape_for_curation_api,
)
from backend.layers.auth.user_info import UserInfo
from backend.layers.business.entities import CollectionMetadataUpdate
from backend.layers.business.exceptions import CollectionUpdateException, InvalidMetadataException
from backend.layers.common.entities import CollectionId, CollectionLinkType, Link
from backend.portal.api.providers import get_business_logic, get_cloudfront_provider


def delete(collection_id: str, token_info: dict, delete_published: bool = False) -> Response:
    user_info = UserInfo(token_info)
    collection_version = get_inferred_collection_version(collection_id)
    is_owner_or_allowed_else_forbidden(collection_version, user_info)
    if collection_version.published_at:
        if user_info.is_cxg_admin() and delete_published:
            get_business_logic().tombstone_collection(CollectionId(collection_id))
            get_cloudfront_provider().create_invalidation_for_index_paths()
        elif not user_info.is_cxg_admin() and delete_published:
            raise ForbiddenHTTPException(
                detail="Only CXG admins can delete a published collection. Roles are granted through Auth0"
            )
        else:
            raise MethodNotAllowedException(detail="Cannot delete a published collection through API.")
    else:
        get_business_logic().delete_collection_version(collection_version)
        get_cloudfront_provider().create_invalidation_for_index_paths()
    return make_response("", 204)


def get(collection_id: str, token_info: dict) -> Response:
    collection_version = get_inferred_collection_version(collection_id)
    user_info = UserInfo(token_info)
    response = reshape_for_curation_api(collection_version, user_info)
    return jsonify(response)


def patch(collection_id: str, body: dict, token_info: dict) -> Response:
    user_info = UserInfo(token_info)

    if "links" in body and not body["links"]:
        raise InvalidParametersHTTPException(detail="If provided, the 'links' array may not be empty")

    collection_version = get_inferred_collection_version(collection_id)
    is_owner_or_allowed_else_forbidden(collection_version, user_info)
    if collection_version.published_at is not None:
        raise MethodNotAllowedException(
            detail="Directly editing a public Collection is not allowed; you must create a revision."
        )

    errors = []
    # Verify DOI
    if (doi_url := body.pop("doi", None)) and (doi_url := doi.curation_get_normalized_doi_url(doi_url, errors)):
        links = body.get("links", [])
        links.append({"link_type": CollectionLinkType.DOI.name, "link_url": doi_url})
        body["links"] = links

    # TODO: dedup
    def _link_from_request(body: dict):
        return Link(
            body.get("link_name"),
            body["link_type"],
            body["link_url"],
        )

    update_links = [_link_from_request(node) for node in body["links"]] if body.get("links") is not None else None

    # Build CollectionMetadataUpdate object
    collection_metadata = CollectionMetadataUpdate(
        body.get("name"),
        body.get("description"),
        body.get("contact_name"),
        body.get("contact_email"),
        update_links,
        body.get("consortia", []),
    )

    # only apply doi updates if update_links is populated
    apply_doi_update = update_links is not None

    # Update the collection
    try:
        get_business_logic().update_collection_version(
            collection_version.version_id, collection_metadata, apply_doi_update=apply_doi_update
        )
    except InvalidMetadataException as ex:
        raise InvalidParametersHTTPException(ext=dict(invalid_parameters=ex.errors)) from None
    except CollectionUpdateException as ex:
        errors.extend(ex.errors)
    if errors:
        raise InvalidParametersHTTPException(ext=dict(invalid_parameters=errors))
    # Make Response
    updated_collection_version = get_business_logic().get_collection_version(collection_version.version_id)

    metadata = asdict(updated_collection_version.metadata)
    metadata.pop("links")

    doi_url, links = extract_doi_from_links(updated_collection_version.metadata.links)

    response = dict(
        **metadata, publisher_metadata=updated_collection_version.publisher_metadata, links=links, doi=doi_url
    )
    return jsonify(response)
