import os
from flask import make_response, jsonify


def get():
    return make_response(jsonify({"Data Portal": os.environ["COMMIT_SHA"]}), 200)
