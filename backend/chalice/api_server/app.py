import os
import sys
from functools import wraps

import connexion
from collections import defaultdict
import chalice
import re
from chalice import Chalice, CORSConfig
import logging

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "chalicelib"))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from corpora.common.utils.db_utils import DbUtils
from corpora.common.utils.s3_utils import generate_file_url
from corpora.common.authorizer import assert_authorized


cors_config = CORSConfig(allow_origin="*", max_age=600, allow_credentials=True)

def requires_auth():
    """
    A decorator for assert_authorized
    :return: 401 or the original function response.
    """

    def decorate(func):
        @wraps(func)
        def call(*args, **kwargs):
            assert_authorized(app.current_request.headers)
            return func(*args, **kwargs)

        return call

    return decorate


def create_flask_app():
    app = connexion.App(f"{os.environ['APP_NAME']}-{os.environ['DEPLOYMENT_STAGE']}")
    swagger_spec_path = os.path.join(pkg_root, "config", f"{os.environ['APP_NAME']}.yml")
    logging.warning(f"Swagger spec path: {swagger_spec_path}")
    app.add_api(swagger_spec_path, validate_responses=False)
    return app.app


def get_chalice_app(flask_app):
    app = Chalice(app_name=flask_app.name)
    flask_app.debug = True
    app.debug = flask_app.debug
    app.log.setLevel(logging.DEBUG)

    def dispatch(*args, **kwargs):
        uri_params = app.current_request.uri_params or {}
        path = app.current_request.context["resourcePath"].format(**uri_params)
        req_body = app.current_request.raw_body if app.current_request._body is not None else None
        with flask_app.test_request_context(
            path=path,
            base_url="https://{}".format(app.current_request.headers["host"]),
            query_string=app.current_request.query_params,
            method=app.current_request.method,
            headers=list(app.current_request.headers.items()),
            data=req_body,
            environ_base=app.current_request.stage_vars,
        ):
            flask_res = flask_app.full_dispatch_request()
        res_headers = dict(flask_res.headers)
        # API Gateway/Cloudfront adds a duplicate Content-Length with a different value (not sure why)
        res_headers.pop("Content-Length", None)
        return chalice.Response(
            status_code=flask_res._status_code,
            headers=res_headers,
            body="".join([c.decode() if isinstance(c, bytes) else c for c in flask_res.response]),
        )

    routes = defaultdict(list)
    for rule in flask_app.url_map.iter_rules():
        routes[re.sub(r"<(.+?)(:.+?)?>", r"{\1}", rule.rule).rstrip("/")] += rule.methods
    for route, methods in routes.items():
        app.route(route, methods=list(set(methods) - {"OPTIONS"}), cors=True)(dispatch)

    with open(os.path.join(pkg_root, "index.html")) as swagger_ui_file_object:
        swagger_ui_html = swagger_ui_file_object.read()

    @app.route("/")
    def serve_swagger_ui():
        return chalice.Response(status_code=200, headers={"Content-Type": "text/html"}, body=swagger_ui_html)

    @app.route("/projects/{project_id}", cors=cors_config)
    def get_project(project_id: str):
        db = DbUtils()
        project = db.query_project(project_id)

        return chalice.Response(
            status_code=200 if project else 404, headers={"Content-Type": "application/json"}, body=project,
        )

    @app.route("/projects/{project_id}/files", cors=cors_config)
    @requires_auth()
    def get_project_files(project_id: str):
        db = DbUtils()
        files = db.query_downloadable_project_files(project_id)
        return chalice.Response(
            status_code=200 if files else 404, headers={"Content-Type": "application/json"}, body=files,
        )

    @app.route("/files/{file_id}", cors=cors_config)
    @requires_auth()
    def get_file(file_id: str):
        db = DbUtils()
        file = db.query_file(file_id)

        download_url = ""
        if file:
            project = db.query_project(file["project_id"])
            file_prefix = f"{project['label']}/matrix.loom"
            download_url = generate_file_url(file_prefix)

        return chalice.Response(
            status_code=200 if file else 404, headers={"Content-Type": "application/json"}, body={"url": download_url},
        )

    return app


app = get_chalice_app(create_flask_app())

