import os
from flask import make_response, jsonify


def get():
    version_info = {"Data Portal": os.getenv("COMMIT_SHA")}
    return make_response(jsonify(version_info), 200)
