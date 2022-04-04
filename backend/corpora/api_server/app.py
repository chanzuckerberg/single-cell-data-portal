from backend.corpora.lambdas.api.v1.authorization import AuthError
import connexion
import json
import logging
import os
from connexion import FlaskApi, ProblemException, problem
from flask import g, jsonify, request
from flask_cors import CORS
from urllib.parse import urlparse

from backend.corpora.common.utils.json import CustomJSONEncoder
from backend.corpora.common.utils.aws import AwsSecret

DEPLOYMENT_STAGE = os.environ["DEPLOYMENT_STAGE"]
APP_NAME = os.environ["APP_NAME"]


def create_flask_app():
    connexion_app = connexion.FlaskApp(f"{APP_NAME}-{DEPLOYMENT_STAGE}", specification_dir="backend/config")
    # From https://github.com/zalando/connexion/issues/346
    connexion_app.app.url_map.strict_slashes = False
    options = {"swagger_ui": True, "swagger_url": "/"}
    dataportal_api = super(connexion.FlaskApp, connexion_app).add_api(
        f"{APP_NAME}.yml", validate_responses=True, options=options
    )
    curator_api = super(connexion.FlaskApp, connexion_app).add_api(
        "curation-api.yml", validate_responses=True, options=options, resolver_error=501
    )
    connexion_app.app.register_blueprint(dataportal_api.blueprint, url_prefix="/", name="data-portal")
    connexion_app.app.register_blueprint(curator_api.blueprint, url_prefix="/curation", name="curation")
    return connexion_app.app


def configure_flask_app(flask_app):
    # configure logging
    gunicorn_logger = logging.getLogger("gunicorn.error")
    flask_app.logger.handlers = gunicorn_logger.handlers
    flask_app.logger.setLevel(gunicorn_logger.level)
    flask_app.debug = False if DEPLOYMENT_STAGE == "prod" else True
    logging.basicConfig(level=gunicorn_logger.level)

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


app = configure_flask_app(create_flask_app())


@app.before_request
def pre_request_logging():
    message = json.dumps(dict(url=request.path, method=request.method, schema=request.scheme))
    app.logger.info(message)


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
