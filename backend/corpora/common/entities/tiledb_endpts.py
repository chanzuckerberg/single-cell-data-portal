# Wraps tiledb_data.py to fit API expected responses
import logging
import time

from flask import make_response, jsonify
from pyrsistent import s
from backend.corpora.common.utils.http_exceptions import TooLargeHTTPException
from backend.corpora.lambdas.api.v1.collection_id.upload import upload_from_link
from backend.corpora.lambdas.api.v1.dataset import post_dataset_asset
from backend.corpora.lambdas.api.v1.collection import normalize_links_with_doi, get_publisher_metadata
from backend.corpora.lambdas.api.v1.authorization import is_user_owner_or_allowed

from backend.corpora.common.utils.http_exceptions import InvalidParametersHTTPException, TooLargeHTTPException, ForbiddenHTTPException, MethodNotAllowedException, NotFoundHTTPException, ServerErrorHTTPException

from backend.corpora.common.entities.tiledb_data import TileDBData

logger = logging.getLogger(__name__)

db = TileDBData()

def create_collection(body: dict, user: str):
    """/v1/collections POST"""
    for link in body.get("links", []):
        link["link_url"] = link["link_url"].strip()

    errors = []
    links, doi = normalize_links_with_doi(body, errors)
    publisher_metadata = {}
    if errors:
        return make_response({"detail": errors}, 400)
    if doi is not None:
        data = get_publisher_metadata(doi)
        publisher_metadata = data if data else {}

    id = db.create_collection(
        dict(
            visibility="PRIVATE",
            name=body["name"].strip(),
            description=body["description"].strip(),
            owner=user,
            links=links,
            contact_name=body["contact_name"].strip(),
            contact_email=body["contact_email"].strip(),
            curator_name=body.get("curator_name", "").strip(),
            publisher_metadata=publisher_metadata
        )
    )
    res = {"collection_id": id}
    return make_response(res, 201)


def list_collections(token_info: dict = None):
    """/v1/collections GET"""
    colls = db.get_all_collections()

    filtered = []
    for coll in colls:
        visibility = coll["visibility"]
        owner = coll["owner"]
        if visibility == "PUBLIC":
            filtered.append(dict(id=coll["id"], created_at=coll["created_at"], visibility=visibility))
        elif is_user_owner_or_allowed(token_info, owner) and visibility != "DELETED":
            filtered.append(
                dict(
                    id=coll["id"],
                    created_at=coll["created_at"],
                    visibility=visibility,
                    revision_of=coll["revision_of"],
                )
            )
    return make_response({"collections": filtered}, 200)


def list_published_collections_compact():
    """/v1/collections/index GET"""
    data = db.get_published_collections()
    colls = [{
        'id': x['id'], 'name': x['name'],
        'publisher_metadata': x['publisher_metadata'],
        'revised_at': x['updated_at']} 
    for x in data]
    return make_response(jsonify(colls))


def delete_collection(collection_id: str, token_info: dict):
    """/v1/collections/{id} DELETE"""
    
    coll = db.get_collection(collection_id)
    if not coll:
        return make_response({"detail": "Collection not found."}, 403)
    elif coll['visibility'] == "DELETED":
        return make_response('', 403)
    if is_user_owner_or_allowed(token_info, coll['owner']) or coll['visibility'] == "PRIVATE":
        db.delete_collection(collection_id)
        return make_response('', 204)
    else:
        return make_response({"detail": "Cannot delete collection"}, 401)


def get_collection(collection_id: str, token_info: dict):
    """/v1/collections/{id} GET"""
    coll = db.get_collection(collection_id)
    if not coll:
        return make_response({"detail": "Collection not found."}, 403)
    elif coll['visibility'] == "DELETED":
        return make_response("", 403)
        
    ids = coll['datasets']
    datasets = db.get_datasets(collection_id, ids)

    coll['datasets'] = datasets
    owner = coll['owner']
    coll["access_type"] = "WRITE" if is_user_owner_or_allowed(token_info, owner) else "READ"
    res = reshape_coll_for_api(coll)
    return make_response(res, 200)


def start_revision(collection_id: str, token_info: dict):
    """/v1/collections/{id} POST"""
    
    coll = db.get_collection(collection_id)
    if not coll:
        return make_response({"detail": "Collection not found."}, 404)
    if len(coll['revision_of']) > 0:
        return make_response({"detail": "Cannot revise a revision."}, 403)
    if is_user_owner_or_allowed(token_info, coll['owner']) and coll['visibility'] == "PRIVATE":
        rev_id = db.create_revision(collection_id)
        res = get_collection(rev_id)
        res["access_type"] = "WRITE"
        return make_response(res, 200)
    else:
        return make_response({"detail": "Cannot revise this collection"}, 403)


def update_collection(collection_id: str, body: dict, token_info: dict):
    """/v1/collections/{id} PUT"""
    coll = db.get_collection(collection_id)
    if not coll:
        return make_response({"detail": "Collection not found."}, 404)
    if not is_user_owner_or_allowed(token_info, coll['owner']):
        make_response({"detail": "Cannot update collection owned by a different user"}, 403)
    if coll['visibility'] != "PRIVATE":
        make_response({"detail": "Cannot update non-private collection."}, 403)

    for key, val in body.items():
        if val != coll[key]:
            if key == "links":
                errors = []
                links, doi = normalize_links_with_doi(body, errors)
                _, old_doi = normalize_links_with_doi(coll, [])
                publisher_metadata = {}
                if errors:
                    return make_response({"detail": errors}, 400)
                if doi is not None and doi != old_doi:
                    data = get_publisher_metadata(doi)
                    publisher_metadata = data if data else {}
                    db.edit_collection(collection_id, "publisher_metadata", publisher_metadata)
                val = links
            db.edit_collection(collection_id, key, val)

    res = get_collection(collection_id, token_info)
    return make_response(res, 200)


def publish_collection(collection_id: str, token_info: dict):
    """/v1/collections/{id}/publish POST"""
    
    coll = db.get_collection(collection_id)
    if not coll:
        return make_response({"detail": "Collection not found."}, 404)
    ids = coll['datasets']
    datasets = db.get_datasets(collection_id, ids)
    if not is_user_owner_or_allowed(token_info, coll['owner']) or len(datasets) == 0:
        return make_response({"detail": "Cannot publish this collection"}, 403)
    db.publish_collection(collection_id)
    res = {
        "collection_id": collection_id,
        "visibility": "PUBLIC"
    }
    return make_response(res, 202)


def upload_dataset(collection_id: str, body: dict, token_info: dict):
    """/v1/collections/{id}/upload-links POST"""
    # get the dataset data from the given url, manage and upload the artifact, etc.
    try:
        dataset_id = upload_from_link(collection_id, token_info, body["url"], curator_tag=body.get("curator_tag"))
        res = {"dataset_id": dataset_id}
        return make_response(res, 202)
    except TooLargeHTTPException:
        return make_response({"detail": "File too large to upload."}, 413)
    except InvalidParametersHTTPException:
        return make_response({"detail": "Invalid parameters."}, 400)
    except ForbiddenHTTPException:
        return make_response({"detail": "Unauthorized to upload to this collection."}, 401)
    except MethodNotAllowedException:
        return make_response({"detail": "Invalid processing status."}, 409)
    except NotFoundHTTPException:
        return make_response({"detail": "Non-existent dataset."}, 404)
    except Exception:
        return make_response({"detail": "Invalid dataset"}, 400)


def replace_dataset(collection_id: str, body: dict, token_info: dict):
    """/v1/collections/{id}/upload-links PUT"""
    try:
        
        dataset_id = body.get("id")
        dataset_id = upload_from_link(
            collection_id,
            token_info,
            body.get("url", body.get("link")),
            dataset_id,
            curator_tag=body.get("curator_tag"),
        )
        db.delete_dataset(collection_id, dataset_id)
        res = {"dataset_id": dataset_id}
        return make_response(res, 202)
    except TooLargeHTTPException:
        return make_response({"detail": "File too large to upload."}, 413)
    except InvalidParametersHTTPException:
        return make_response({"detail": "Invalid parameters."}, 400)
    except ForbiddenHTTPException:
        return make_response({"detail": "Unauthorized to upload to this collection."}, 401)
    except MethodNotAllowedException:
        return make_response({"detail": "Invalid processing status."}, 409)
    except NotFoundHTTPException:
        return make_response({"detail": "Non-existent dataset."}, 404)


def list_published_datasets():
    """/v1/datasets/index GET"""
    
    datasets = db.get_published_datasets()
    return make_response(jsonify(datasets), 200)


def get_dataset_metadata(explorer_url: str):
    """/v1/datasets/meta GET"""
    # TODO
    return


def delete_dataset(collection_id: str, dataset_id: str, token_info: dict):
    """/v1/datasets/{id} DELETE"""
    
    coll = db.get_collection(collection_id)
    if not coll:
        return make_response({"detail": "Collection not found."}, 404)
    if not is_user_owner_or_allowed(token_info, coll['owner']):
        return make_response({"detail": "Unauthorized to delete this dataset."}, 403)
    if coll['visibility'] == "PUBLIC":
        return make_response({"detail": "Cannot delete a public dataset."}, 405)
    if dataset_id not in coll['datasets']:
        return make_response({"detail": "Dataset not found."}, 405)
    db.delete_dataset(collection_id, dataset_id)
    return make_response('', 202)


def request_asset_download(dataset_id: str, asset_id: str):
    """/v1/datasets/{dataset_id}/asset/{asset_id} POST"""
    try:
        res = post_dataset_asset(dataset_id, asset_id)
        return make_response(res, 200)
    except NotFoundHTTPException:
        return make_response({"detail": "Dataset or asset not found."}, 404)
    except ServerErrorHTTPException:
        return make_response({"detail": "Could not retrieve data or generate URL."}, 404)


def list_dataset_assets(dataset_id: str):
    """/v1/datasets/{dataset_id}/assets GET"""
    
    dataset = db.get_dataset(dataset_id)
    if not dataset:
        return make_response({"detail": "Dataset not found."}, 404)
    assets = dataset["dataset_assets"]
    return make_response(jsonify(assets), 200)


def get_dataset_upload_status(collection_id: str, dataset_id: str, token_info: dict):
    """/v1/datasets/{dataset_id}/status GET"""
    
    coll = db.get_collection(collection_id)
    if not coll:
        return make_response({"detail": "Collection not found."}, 404)
    if not is_user_owner_or_allowed(token_info, coll['owner']):
        return make_response({"detail": "Unauthorized to read this dataset."}, 403)
    dataset = db.get_dataset(dataset_id)
    if not dataset:
        return make_response({"detail": "Dataset not found."}, 404)
    status = dataset['processing_status']
    return make_response(status, 200)


def reshape_coll_for_api(result: dict):
    """Modify some parts of a collection object to fit what the API and frontend expect. """
    # TODO: Ideally, we wouldn't need this kind of adaptor at all
    for hidden in ["user", "owner"]:
        result.pop(hidden, None)
    result["links"] = [
        dict(link_url=link["link_url"], link_name=link.get("link_name", ""), link_type=link["link_type"])
        for link in result["links"]
    ]
    datasets = []
    for dataset in result["datasets"]:
        dataset["dataset_deployments"] = []
        explorer_url = dataset.pop("explorer_url", None)
        if explorer_url:
            dataset["dataset_deployments"].append({"url": explorer_url})
        datasets.append(dataset)

        from backend.corpora.common.entities import Dataset
        Dataset.transform_sex_for_schema_2_0_0(dataset)
        # Dataset.transform_organism_for_schema_2_0_0(dataset)

    result["datasets"] = datasets
    return result