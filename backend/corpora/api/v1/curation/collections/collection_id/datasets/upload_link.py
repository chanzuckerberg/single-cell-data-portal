from flask import make_response

from backend.corpora.api.v1.collection_id.upload import upload_from_link


def put(collection_id: str, body: dict, token_info: dict):
    dataset_id = upload_from_link(
        collection_id,
        token_info,
        body.get("url", body.get("link")),
        body.get("id"),
        curator_tag=body.get("curator_tag"),
    )
    return make_response({"dataset_id": dataset_id}, 202)
