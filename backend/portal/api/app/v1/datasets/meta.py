from flask import g, make_response, jsonify

from backend.api_server.db import dbconnect
from backend.common.corpora_orm import DatasetArtifactFileType
from backend.common.entities import Dataset
from backend.common.utils.http_exceptions import NotFoundHTTPException


@dbconnect
def get(url: str):
    db_session = g.db_session
    dataset = Dataset.get_by_explorer_url(db_session, url)
    if not dataset:
        raise NotFoundHTTPException()
    artifact = dataset.get_most_recent_artifact(filetype=DatasetArtifactFileType.CXG)
    s3_uri = artifact.s3_uri if artifact else None

    dataset_identifiers = {
        "s3_uri": s3_uri,
        "dataset_id": dataset.id,
        "collection_id": dataset.collection_id,
        "collection_visibility": dataset.collection.visibility,
        "tombstoned": dataset.tombstone,
    }
    return make_response(jsonify(dataset_identifiers), 200)
