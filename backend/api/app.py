import json
import logging
import os
from typing import Dict
from urllib.parse import urlparse

import connexion
from connexion import FlaskApi, ProblemException, problem
from flask import g, jsonify
from flask_cors import CORS
from swagger_ui_bundle import swagger_ui_3_path

from backend.api.data_portal.common.utils.aws import AwsSecret
from backend.api.data_portal.common.utils.json import CustomJSONEncoder
from backend.api.data_portal.lambdas.api.v1.authorization import AuthError

DEPLOYMENT_STAGE = os.environ["DEPLOYMENT_STAGE"]
APP_NAME = os.environ["APP_NAME"]


def create_flask_app(apis: Dict):
    connexion_app = connexion.FlaskApp(
        f"{APP_NAME}-{DEPLOYMENT_STAGE}", specification_dir=f"{os.path.dirname(__file__)}"
    )
    # From https://github.com/zalando/connexion/issues/346
    connexion_app.app.url_map.strict_slashes = False

    # Add each API under its own base path
    for base_path, spec_file in apis.items():
        connexion_app.add_api(
            spec_file,
            validate_responses=True,
            base_path=f"/{base_path}",
            options={
                "serve_spec": True,
                "swagger_path": swagger_ui_3_path,
                "swagger_ui": True,
                "swagger_url": None,
                "verbose": True,
            },
        )

    return connexion_app.app


def configure_flask_app(flask_app):
    # configure logging
    gunicorn_logger = logging.getLogger("gunicorn.error")
    flask_app.logger.handlers = gunicorn_logger.handlers
    flask_app.logger.setLevel(gunicorn_logger.level)
    flask_app.debug = False if DEPLOYMENT_STAGE == "prod" else True

    # set the flask secret key, needed for session cookies
    flask_secret_key = "OpenSesame"
    allowed_origins = []
    deployment_stage = os.environ["DEPLOYMENT_STAGE"]
    if deployment_stage not in ["prod"]:
        allowed_origins.extend([r"http://.*\.corporanet\.local:\d+", r"^http://localhost:\d+"])
    # TODO: Move the FRONTEND_URL env var case to AppConfigPropertiesSource as EnvConfigPropertiesSource.
    # Is this even used?
    if os.getenv("FRONTEND_URL"):
        allowed_origins.append(os.getenv("FRONTEND_URL"))
    if deployment_stage != "test":  # pragma: no cover
        secret_name = f"corpora/backend/{deployment_stage}/auth0-secret"
        auth_secret = json.loads(AwsSecret(secret_name).value)
        if auth_secret:
            flask_secret_key = auth_secret.get("flask_secret_key", flask_secret_key)
            frontend = auth_secret.get("redirect_to_frontend", None)
            if frontend:
                if frontend.endswith("/"):
                    frontend = frontend[:-1]
                frontend_parse = urlparse(frontend)
                allowed_origins.append(f"{frontend_parse.scheme}://{frontend_parse.netloc}")
    flask_app.logger.info(f"CORS allowed_origins: {allowed_origins}")
    CORS(flask_app, max_age=600, supports_credentials=True, origins=allowed_origins, allow_headers=["Content-Type"])

    # FIXME, enforce that the flask_secret_key is found once all secrets are setup for all environments
    require_secure_cookies = not bool(os.getenv("DEV_MODE_COOKIES"))
    flask_app.config.update(
        SECRET_KEY=flask_secret_key,
        SESSION_COOKIE_SECURE=require_secure_cookies,
        SESSION_COOKIE_HTTPONLY=True,
        SESSION_COOKIE_SAMESITE="Lax",
    )
    flask_app.json_encoder = CustomJSONEncoder
    return flask_app


class InterceptRequestMiddleware:
    def __init__(self, wsgi_app):
        self.wsgi_app = wsgi_app

    def __call__(self, environ, start_response):
        environ["HTTP_CXGPUBLIC"] = "dummy"  # TODO: What is this? Is it used?
        return self.wsgi_app(environ, start_response)


# TODO: Make this config automatic, and document how to add a new API
apis = {"wmg": "wmg/wmg-api.yml", "dp": "data_portal/data-portal-api.yml"}
app = configure_flask_app(create_flask_app(apis))
app.wsgi_app = InterceptRequestMiddleware(app.wsgi_app)


@app.route("/")
def apis_landing_page():
    # TODO: use jinja2 template to render this
    links = [f'<a href="{api_name}/ui/">{api_name}</a></br>' for api_name in apis.keys()]
    return f"""
    <html>
      <head><title>cellxgene Platform APIs</title></head>
      <body><h1>cellxgene Platform APIs</h1>{"".join(links)}</body>
    </html>
    """


@app.teardown_appcontext
def close_db(e=None):
    g.pop("db_session", None)


@app.errorhandler(AuthError)
def handle_auth_error(ex):
    response = jsonify(ex.error)
    response.status_code = ex.status_code
    return response


@app.errorhandler(ProblemException)
def handle_corpora_error(exception):
    return FlaskApi.get_response(
        problem(
            exception.status,
            exception.title,
            exception.detail,
            exception.type,
            exception.instance,
            exception.headers,
            exception.ext,
        )
    )


if __name__ == "__main__":
    app.run(host="0.0.0.0", debug=True)
