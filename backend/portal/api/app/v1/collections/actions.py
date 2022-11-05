from typing import Optional

from flask import make_response, jsonify, g

from backend.portal.api.app.v1.common import portal_get_normalized_doi_url
from backend.portal.api.collections_common import create_collection_common
from backend.common.corpora_orm import DbCollection, CollectionVisibility, ProjectLinkType
from backend.common.entities import Collection
from backend.portal.api.app.v1.authorization import is_user_owner_or_allowed
from backend.api_server.db import dbconnect


@dbconnect
def get(from_date: int = None, to_date: int = None, token_info: Optional[dict] = None):
    db_session = g.db_session
    all_collections = Collection.list_attributes_in_time_range(
        db_session,
        from_date=from_date,
        to_date=to_date,
        list_attributes=[
            DbCollection.id,
            DbCollection.visibility,
            DbCollection.owner,
            DbCollection.created_at,
            DbCollection.revision_of,
        ],
    )

    collections = []
    for coll_dict in all_collections:
        visibility = coll_dict["visibility"]
        owner = coll_dict["owner"]
        if visibility == CollectionVisibility.PUBLIC:
            collections.append(dict(id=coll_dict["id"], created_at=coll_dict["created_at"], visibility=visibility.name))
        elif is_user_owner_or_allowed(token_info, owner):
            collections.append(
                dict(
                    id=coll_dict["id"],
                    created_at=coll_dict["created_at"],
                    visibility=visibility.name,
                    revision_of=coll_dict["revision_of"],
                )
            )

    result = {"collections": collections}
    if from_date:
        result["from_date"] = from_date
    if to_date:
        result["to_date"] = to_date

    return make_response(jsonify(result), 200)


@dbconnect
def get_collections_index():
    # TODO (ebezzi): this is very similar to `get` above. Eventually they should be consolidated
    db_session = g.db_session

    filtered_collection = Collection.list_attributes_in_time_range(
        db_session,
        filters=[DbCollection.visibility == CollectionVisibility.PUBLIC],
        list_attributes=[
            DbCollection.id,
            DbCollection.name,
            DbCollection.published_at,
            DbCollection.revised_at,
            DbCollection.publisher_metadata,
        ],
    )

    # Remove entries where the value is None
    updated_collection = []
    for d in filtered_collection:
        updated_collection.append({k: v for k, v in d.items() if v is not None})

    return make_response(jsonify(updated_collection), 200)


def get_doi_link_node(body: dict, errors: list) -> Optional[dict]:
    links = body.get("links", [])
    dois = [link for link in links if link["link_type"] == ProjectLinkType.DOI.name]

    if not dois:
        return None

    # Verify that a single DOI exists
    if len(dois) > 1:
        errors.append({"link_type": ProjectLinkType.DOI.name, "reason": "Can only specify a single DOI"})
        return None

    return dois[0]


def post(body: dict, user: str):
    errors = []
    doi_url = None
    if doi_node := get_doi_link_node(body, errors):
        if doi_url := portal_get_normalized_doi_url(doi_node, errors):
            doi_node["link_url"] = doi_url
    collection_id = create_collection_common(body, user, doi_url, errors)
    return make_response(jsonify({"collection_id": collection_id}), 201)
