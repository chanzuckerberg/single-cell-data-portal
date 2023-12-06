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

    cell_onthology_id = body["cell_onthology_id"]
    references = body["references"]
    description = body["description"]

    if not re.match(r"^CL_\d{7}$", cell_onthology_id):
        raise ForbiddenHTTPException("Invalid cell_onthology_id format. Example: CL_0000030")

    for reference in references:
        try:
            requests.get(reference)
        except Exception as e:
            raise ForbiddenHTTPException(f"Invalid url {reference}") from e

    file_content = {"description": description, "references": references}

    if os.getenv("DEPLOYMENT_STAGE") == "rdev" or os.getenv("DEPLOYMENT_STAGE") == "test":
        env = "dev"
    else:
        env = os.getenv("DEPLOYMENT_STAGE")

    file_name = f"{cell_onthology_id}.json"
    key_name = f"validated_descriptions/{file_name}"
    bucket_name = f"cellguide-data-public-{env}"

    s3_provider = S3Provider()
    file_content_str = json.dumps(file_content)
    s3_provider.put_object(bucket_name, key_name, file_content_str)

    file_content["cell_onthology_id"] = cell_onthology_id
    return make_response(jsonify(file_content), 201)
