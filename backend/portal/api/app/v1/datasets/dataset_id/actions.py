from flask import g

from backend.common.entities import Dataset
from backend.api_server.db import dbconnect
from backend.portal.api.collections_common import delete_dataset_common


@dbconnect
def delete(dataset_id: str, token_info: dict):
    """
    Deletes an existing dataset or cancels an in progress upload.
    """
    db_session = g.db_session
    dataset = Dataset.get(db_session, dataset_id, include_tombstones=True)
    delete_dataset_common(db_session, dataset, token_info)
    return "", 202
