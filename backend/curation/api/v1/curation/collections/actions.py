from flask import jsonify, make_response

from backend.common.utils.http_exceptions import ForbiddenHTTPException, InvalidParametersHTTPException
from backend.curation.api.v1.curation.collections.common import reshape_for_curation_api
from backend.layers.auth.user_info import UserInfo
from backend.layers.business.entities import CollectionQueryFilter
from backend.layers.business.exceptions import CollectionCreationException, InvalidMetadataException
from backend.layers.common import doi
from backend.layers.common.entities import CollectionLinkType, CollectionMetadata, CollectionVisibility, Link
from backend.portal.api.providers import get_business_logic


def get(visibility: str, token_info: dict, curator: str = None):  # type: ignore
    """
    Collections index endpoint for Curation API. Only return Collection data for which the curator is authorized.
    :param visibility: the CollectionVisibility in string form
    :param token_info: access token info
    :param curator: the name of the curator to filter the collections by.

    @return: Response
    """
    user_info = UserInfo(token_info)
    filters = {}
    if visibility == CollectionVisibility.PRIVATE.name:
        filters["is_published"] = False
        if user_info.is_none():
            raise ForbiddenHTTPException(detail="Not authorized to query for PRIVATE collection.")
        else:
            if not user_info.is_super_curator():  # A super curator and don't need to filter by owner.
                filters["owner"] = user_info.user_id
    else:
        filters["is_published"] = True

    if curator:
        if not user_info.is_super_curator():
            raise ForbiddenHTTPException(detail="Not authorized to use the curator query parameter.")
        else:
            filters["curator_name"] = curator  # type: ignore

    print(filters)

    resp_collections = []
    for collection_version in get_business_logic().get_collections(CollectionQueryFilter(**filters)):  # type: ignore
        resp_collection = reshape_for_curation_api(collection_version, user_info, preview=True)
        resp_collections.append(resp_collection)
    return jsonify(resp_collections)


def post(body: dict, token_info: dict):
    # Extract DOI into link
    errors = []  # type: ignore
    if (doi_url := body.get("doi")) and (doi_url := doi.curation_get_normalized_doi_url(doi_url, errors)):
        links = body.get("links", [])
        links.append({"link_type": CollectionLinkType.DOI.name, "link_url": doi_url})
        body["links"] = links

    # Build CollectionMetadata object
    links = [Link(link.get("link_name"), link["link_type"], link["link_url"]) for link in body.get("links", [])]
    metadata = CollectionMetadata(
        body["name"],
        body["description"],
        body["contact_name"],
        body["contact_email"],
        links,
        body.get("consortia", []),
    )

    try:
        version = get_business_logic().create_collection(token_info["sub"], token_info["curator_name"], metadata)
    except InvalidMetadataException as ex:
        errors.extend(ex.errors)  # type: ignore
    except CollectionCreationException as ex:
        errors.extend(ex.errors)  # type: ignore
    if errors:
        raise InvalidParametersHTTPException(ext=dict(invalid_parameters=errors))
    return make_response(jsonify({"collection_id": version.collection_id.id}), 201)
