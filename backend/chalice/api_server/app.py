import logging
import os
import re
import sys
from collections import defaultdict
from functools import wraps

import chalice
import connexion
from chalice import Chalice, CORSConfig

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "chalicelib"))  # noqa
sys.path.insert(0, pkg_root)  # noqa

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
    app = connexion.FlaskApp(f"{os.environ['APP_NAME']}-{os.environ['DEPLOYMENT_STAGE']}")
    swagger_spec_path = os.path.join(pkg_root, "config", f"{os.environ['APP_NAME']}.yml")
    app.add_api(swagger_spec_path, validate_responses=True)
    return app.app


def get_chalice_app(flask_app):
    app = Chalice(app_name=flask_app.name)
    flask_app.debug = True
    app.debug = flask_app.debug
    app.log.setLevel(logging.DEBUG)

    def clean_entry_for_logging(entry):
        log = entry.to_dict()
        log.pop("body")
        return log

    def dispatch(*args, **kwargs):
        app.log.info(f"Request: {clean_entry_for_logging(app.current_request)}")

        uri_params = app.current_request.uri_params or {}
        resource_path = app.current_request.context["resourcePath"].format(**uri_params)
        req_body = app.current_request.raw_body if app.current_request._body is not None else None

        with flask_app.test_request_context(
            path=resource_path,
            base_url="https://{}".format(app.current_request.headers["host"]),
            query_string=app.current_request.query_params,
            method=app.current_request.method,
            headers=list(app.current_request.headers.items()),
            data=req_body,
            environ_base=app.current_request.stage_vars,
        ):
            flask_res = flask_app.full_dispatch_request()

        response_headers = dict(flask_res.headers)
        response_headers.update({"X-AWS-REQUEST-ID": app.lambda_context.aws_request_id})

        chalice_response = chalice.Response(
            status_code=flask_res._status_code,
            headers=response_headers,
            body="".join([c.decode() if isinstance(c, bytes) else c for c in flask_res.response]),
        )

        app.log.info(f"Response: {clean_entry_for_logging(chalice_response)}")

        return chalice_response

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

    return app


app = get_chalice_app(create_flask_app())
