from flask import make_response

from .....common.utils.db_utils import db_session_manager


def delete(collection_uuid: str):
    with db_session_manager() as session:
        return make_response({}, 202)
