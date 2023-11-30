import json
import os
import re

import requests
from flask import jsonify, make_response

from backend.cellguide.pipeline.providers.s3_provider import S3Provider
from backend.common.utils.http_exceptions import ForbiddenHTTPException
from backend.layers.auth.user_info import UserInfo


def post(body: dict, token_info: dict):
    try:
        UserInfo(token_info).is_cxg_admin()
    except Exception as e:
        raise ForbiddenHTTPException("Unauthorized") from e

    cell_onthology_id = body["cell_onthology_id"]
    references = body["references"]
    description = body["description"]

    if not re.match(r"^CL_\d{7}$", cell_onthology_id):
        raise ForbiddenHTTPException("Invalid cell_onthology_id format. Example: CL_0000030")

    for reference in references:
        try:
            requests.get(reference)
        except:
            raise ForbiddenHTTPException("Invalid url")

    file_content = {"description": description, "references": references}

    env = os.getenv("DEPLOYMENT_STAGE")
    file_name = f"{cell_onthology_id}.json"
    key_name = f"validated_descriptions/{file_name}"
    bucket_name = f"cellguide-data-public-{env}"
    with open(file_name, "w") as f:
        json.dump(file_content, f)

    s3_provider = S3Provider()
    s3_provider.upload_file(file_name, bucket_name, key_name, {})

    os.remove(file_name)

    file_content["cell_onthology_id"] = cell_onthology_id
    return make_response(jsonify(file_content), 201)
