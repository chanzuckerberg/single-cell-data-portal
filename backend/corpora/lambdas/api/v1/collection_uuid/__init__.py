from flask import make_response


def delete(collection_uuid: str):
    return make_response({}, 202)
