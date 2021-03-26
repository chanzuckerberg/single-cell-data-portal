from collections import defaultdict

import chalice
import connexion
import json
import logging
import os
import re
import sys
from chalice import Chalice
from connexion import FlaskApi, ProblemException, problem
from flask import g
from flask_cors import CORS
from urllib.parse import urlparse

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "chalicelib"))  # noqa
sys.path.insert(0, pkg_root)  # noqa


from corpora.common.utils.json import CustomJSONEncoder
from corpora.common.utils.aws import AwsSecret
from corpora.common.utils.db_session import db_session_manager


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

    # set the flask secret key, needed for session cookies
    flask_secret_key = "OpenSesame"
    deployment = os.environ["DEPLOYMENT_STAGE"]
    allowed_origins = []
    if deployment not in ["prod"]:
        allowed_origins.extend([r"http://.*\.corporanet\.local:\d+", r"^http://localhost:\d+"])
    if os.getenv("FRONTEND_URL"):
        allowed_origins.append(os.getenv("FRONTEND_URL"))
    if deployment != "test":  # pragma: no cover
        secret_name = f"corpora/backend/{os.environ['DEPLOYMENT_STAGE']}/auth0-secret"
        auth_secret = json.loads(AwsSecret(secret_name).value)
        if auth_secret:
            flask_secret_key = auth_secret.get("flask_secret_key", flask_secret_key)
            frontend = auth_secret.get("redirect_to_frontend", None)
            if frontend:
                if frontend.endswith("/"):
                    frontend = frontend[:-1]
                frontend_parse = urlparse(frontend)
                allowed_origins.append(f"{frontend_parse.scheme}://{frontend_parse.netloc}")
    CORS(flask_app, max_age=600, supports_credentials=True, origins=allowed_origins, allow_headers=["Content-Type"])
    app.log.info(f"CORS allowed_origins: {allowed_origins}")

    # FIXME, enforce that the flask_secret_key is found once all secrets are setup for all environments
    require_secure_cookies = True
    if os.getenv("DEV_MODE_COOKIES"):
        require_secure_cookies = False
    flask_app.config.update(
        SECRET_KEY=flask_secret_key,
        SESSION_COOKIE_SECURE=require_secure_cookies,
        SESSION_COOKIE_HTTPONLY=True,
        SESSION_COOKIE_SAMESITE="Lax",
    )

    def clean_entry_for_logging(entry):
        log = entry.to_dict()
        log.pop("body", None)
        return log

    def dispatch(*args, **kwargs):
        app.log.info(f"Request: {clean_entry_for_logging(app.current_request)}")
        uri_params = app.current_request.uri_params or {}
        resource_path = app.current_request.context["resourcePath"].format(**uri_params)
        req_body = app.current_request.raw_body if app.current_request._body is not None else None

        # Must convert the chalice.MultiDict into a list of tuples. Chalice returns chalice.Multidict which is
        # incompatible with the werkzeug.MultiDict expected by Flask.
        query_string = list(app.current_request.query_params.items()) if app.current_request.query_params else None

        # set dummy auth token value for optional security endpoints
        headers = [*app.current_request.headers.items(), ("cxgpublic", "dummy")]
        host = app.current_request.headers.get("host")
        with flask_app.test_request_context(
            path=resource_path,
            base_url="https://{}".format(host) if host else None,
            query_string=query_string,
            method=app.current_request.method,
            headers=headers,
            data=req_body,
            environ_base=app.current_request.stage_vars,
        ):
            with db_session_manager() as db_session:
                g.db_session = db_session
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
        app.route(route, methods=list(set([*methods, "OPTIONS"])))(dispatch)

    with open(os.path.join(pkg_root, "index.html")) as swagger_ui_file_object:
        swagger_ui_html = swagger_ui_file_object.read()

    @app.route("/", methods=["GET", "HEAD"])
    def serve_swagger_ui():
        return chalice.Response(
            status_code=200,
            headers={"Content-Type": "text/html", "X-AWS-REQUEST-ID": app.lambda_context.aws_request_id},
            body=swagger_ui_html,
        )

    flask_app.json_encoder = CustomJSONEncoder

    @flask_app.errorhandler(ProblemException)
    def handle_corpora_error(exception):
        response = problem(
            exception.status,
            exception.title,
            exception.detail,
            exception.type,
            exception.instance,
            exception.headers,
            exception.ext,
        )
        response.headers["X-AWS-REQUEST-ID"] = app.lambda_context.aws_request_id
        return FlaskApi.get_response(response)

    @flask_app.teardown_appcontext
    def close_db(e=None):
        g.pop("db_session", None)

    return app


app = get_chalice_app(create_flask_app())
