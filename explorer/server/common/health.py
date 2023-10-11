from http import HTTPStatus

from flask import jsonify, make_response

from server.common.utils.data_locator import DataLocator
from server.version import __version__ as cellxgene_version


def _is_accessible(path, config):
    try:
        dl = DataLocator(path, region_name=config.server__data_locator__s3_region_name)
        return dl.exists()
    except RuntimeError:
        return False


def health_check(config):
    """
    simple health check - return HTTP response.
    See https://tools.ietf.org/id/draft-inadarei-api-health-check-01.html
    """
    health = {"status": None, "version": "1", "releaseID": cellxgene_version}

    dataroot_paths: list[str] = [
        dataroot_value["dataroot"] for dataroot_value in config.server__multi_dataset__dataroots.values()
    ]
    checks: bool = all([_is_accessible(dataroot_path, config) for dataroot_path in dataroot_paths])

    health["status"] = "pass" if checks else "fail"
    code = HTTPStatus.OK if health["status"] == "pass" else HTTPStatus.BAD_REQUEST
    response = make_response(jsonify(health), code)
    response.headers["Content-Type"] = "application/health+json"
    return response
