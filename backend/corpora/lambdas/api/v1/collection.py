import sqlalchemy
from typing import Optional

from flask import make_response, jsonify, g

from ....common.corpora_orm import DbCollection, CollectionVisibility
from ....common.entities import Collection
from ....common.utils.exceptions import ForbiddenHTTPException, ConflictException
from ....api_server.db import dbconnect


@dbconnect
def get_collections_list(from_date: int = None, to_date: int = None, user: Optional[str] = None):
    db_session = g.db_session
    all_collections = Collection.list_attributes_in_time_range(
        db_session,
        from_date=from_date,
        to_date=to_date,
        list_attributes=[DbCollection.id, DbCollection.visibility, DbCollection.owner, DbCollection.created_at],
    )

    collections = []
    for coll_dict in all_collections:
        visibility = coll_dict["visibility"]
        owner = coll_dict["owner"]
        if visibility == CollectionVisibility.PUBLIC or (user and user == owner):
            collections.append(dict(id=coll_dict["id"], created_at=coll_dict["created_at"], visibility=visibility.name))

    result = {"collections": collections}
    if from_date:
        result["from_date"] = from_date
    if to_date:
        result["to_date"] = to_date

    return make_response(jsonify(result), 200)


@dbconnect
def get_collection_details(collection_uuid: str, visibility: str, user: str):
    db_session = g.db_session
    collection = Collection.get_collection(db_session, collection_uuid, visibility)
    if not collection:
        raise ForbiddenHTTPException()
    get_tombstone_datasets = user == collection.owner and collection.visibility == CollectionVisibility.PRIVATE
    result = collection.reshape_for_api(get_tombstone_datasets)
    result["access_type"] = "WRITE" if user == collection.owner else "READ"
    return make_response(jsonify(result), 200)


@dbconnect
def post_collection_revision(collection_uuid: str, user: str):
    db_session = g.db_session
    collection = Collection.get_collection(db_session, collection_uuid, CollectionVisibility.PUBLIC.name, owner=user)
    if not collection:
        raise ForbiddenHTTPException()
    try:
        collection_revision = collection.revision()
    except sqlalchemy.exc.IntegrityError as ex:
        db_session.rollback()
        raise ConflictException() from ex
    result = collection_revision.reshape_for_api()

    result["access_type"] = "WRITE"
    return make_response(jsonify(result), 201)


@dbconnect
def create_collection(body: object, user: str):
    db_session = g.db_session
    collection = Collection.create(
        db_session,
        visibility=CollectionVisibility.PRIVATE,
        name=body["name"],
        description=body["description"],
        owner=user,
        links=body["links"],
        contact_name=body["contact_name"],
        contact_email=body["contact_email"],
        data_submission_policy_version=body["data_submission_policy_version"],
    )

    return make_response(jsonify({"collection_uuid": collection.id}), 201)


def get_collection_dataset(dataset_uuid: str):
    raise NotImplementedError



@dbconnect
def delete_collection(collection_uuid: str, user: str):
    db_session = g.db_session
    priv_collection = Collection.get_collection(
        db_session, collection_uuid, CollectionVisibility.PRIVATE.name, owner=user, include_tombstones=True
    )
    if priv_collection:
        if not priv_collection.tombstone:
            priv_collection.delete()
        return "", 204
    else:
        return "", 403


@dbconnect
def update_collection(collection_uuid: str, body: dict, user: str):
    db_session = g.db_session
    collection = Collection.get_collection(db_session, collection_uuid, CollectionVisibility.PRIVATE.name, owner=user)
    if not collection:
        raise ForbiddenHTTPException()
    collection.update(**body)
    result = collection.reshape_for_api(tombstoned_datasets=True)
    result["access_type"] = "WRITE"
    return make_response(jsonify(result), 200)
