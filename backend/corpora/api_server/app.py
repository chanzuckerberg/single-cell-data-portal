import json
import logging
import os
from urllib.parse import urlparse

import connexion
from connexion import FlaskApi, ProblemException, problem
from flask import g, request
from flask_cors import CORS
from swagger_ui_bundle import swagger_ui_path

from backend.corpora.common.utils.aws import AwsSecret
from backend.corpora.common.utils.json import CustomJSONEncoder
from backend.corpora.lambdas.api.v1.authorization import AuthError

DEPLOYMENT_STAGE = os.environ["DEPLOYMENT_STAGE"]
APP_NAME = os.environ["APP_NAME"]


def create_flask_app():
    connexion_app = connexion.FlaskApp(f"{APP_NAME}-{DEPLOYMENT_STAGE}", specification_dir="backend/config")

    # From https://github.com/zalando/connexion/issues/346
    connexion_app.app.url_map.strict_slashes = False

    def add_api(base_path, spec_file):
        api_base_paths.append(base_path)
        connexion_app.add_api(
            spec_file,
            validate_responses=True,
            base_path=f"/{base_path}",
            resolver_error=501,
            options={
                "serve_spec": True,
                "swagger_path": swagger_ui_path,
                "swagger_ui": True,
                "swagger_url": None,
                "verbose": True,
            },
        )

    add_api(base_path="/dp", spec_file="corpora-api.yml")
    add_api(base_path="/curation", spec_file="curation-api.yml")
    add_api(base_path="/wmg", spec_file="wmg-api.yml")

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


api_base_paths = []

app = configure_flask_app(create_flask_app())


@app.route("/")
def apis_landing_page() -> str:
    """
    Render a page that displays links to all APIs
    """
    # TODO: use jinja2 template to render this
    links = [f'<a href="{base_path}/ui/">{base_path}</a></br>' for base_path in api_base_paths]
    return f"""
    <html>
      <head><title>cellxgene Platform APIs</title></head>
      <body><h1>cellxgene Platform APIs</h1>{"".join(links)}</body>
    </html>
    """


@app.before_request
def pre_request_logging():
    message = json.dumps(dict(url=request.path, method=request.method, schema=request.scheme))
    app.logger.info(message)


@app.teardown_appcontext
def close_db(e=None):
    g.pop("db_session", None)


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
