from flask import g

from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.http_exceptions import InvalidParametersHTTPException
from backend.corpora.lambdas.api.v1.dataset import delete_dataset_common


@dbconnect
def delete_dataset(token_info: dict, collection_uuid: str, curator_tag: str = None, dataset_uuid = None):
    db_session = g.db_session
    if dataset_uuid:
        dataset = Dataset.get(db_session, dataset_uuid, include_tombstones=True)
    elif curator_tag:
        dataset = Dataset.get_dataset_from_curator_tag(db_session, collection_uuid, curator_tag)
    else:
        raise InvalidParametersHTTPException()
    delete_dataset_common(db_session, dataset, token_info)
    return "", 202
