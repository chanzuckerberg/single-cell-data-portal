import json
import os
import re

import requests
from flask import jsonify, make_response

from backend.common.utils.http_exceptions import ForbiddenHTTPException
from backend.layers.auth.user_info import UserInfo
from backend.layers.thirdparty.s3_provider import S3Provider


def post(body: dict, token_info: dict):
    try:
        UserInfo(token_info).is_cxg_admin()
    except Exception as e:
        raise ForbiddenHTTPException("Unauthorized") from e

    cell_ontology_term_id = body["cell_ontology_term_id"]
    references = body["references"]
    description = body["description"]

    if not re.match(r"^CL_\d{7}$", cell_ontology_term_id):
        raise ForbiddenHTTPException("Invalid cell_ontology_term_id format. Example: CL_0000030")

    for reference in references:
        validate_url(reference)

    is_rdev = os.getenv("DEPLOYMENT_STAGE") == "rdev"
    stack_name = os.getenv("REMOTE_DEV_PREFIX")
    rdev_suffix = f"env-rdev-cellguide/{stack_name}/" if is_rdev else ""
    env = "dev" if os.getenv("DEPLOYMENT_STAGE") in ["rdev", "test"] else os.getenv("DEPLOYMENT_STAGE")

    file_name = f"{cell_ontology_term_id}.json"
    key_name = f"{rdev_suffix}validated_descriptions/{file_name}"
    bucket_name = f"cellguide-data-public-{env}"
    file_content = {"description": description, "references": references}
    file_content_str = json.dumps(file_content)

    s3_provider = S3Provider()
    s3_provider.put_object(bucket_name, key_name, file_content_str)

    file_content["cell_ontology_term_id"] = cell_ontology_term_id
    return make_response(jsonify(file_content), 201)


def validate_url(reference):
    try:
        requests.get(reference)
    except Exception as e:
        raise ForbiddenHTTPException(f"Invalid url {reference}") from e
