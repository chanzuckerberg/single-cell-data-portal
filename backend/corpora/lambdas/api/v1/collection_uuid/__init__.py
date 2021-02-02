from flask import make_response

from .....common.utils.db_utils import db_session


@db_session()
def delete(collection_uuid: str):
    return make_response({}, 202)
