import connexion
import json
import logging
import os
from connexion import FlaskApi, ProblemException, problem
from flask import g, Response
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
    swagger_spec_path = f"{APP_NAME}.yml"
    connexion_app.add_api(swagger_spec_path, validate_responses=True)
    return connexion_app.app


def configure_flask_app(flask_app):
    # configure logging
    gunicorn_logger = logging.getLogger('gunicorn.error')
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


class InterceptRequestMiddleware:
    def __init__(self, wsgi_app):
        self.wsgi_app = wsgi_app

    def __call__(self, environ, start_response):
        environ["HTTP_CXGPUBLIC"] = "dummy"
        return self.wsgi_app(environ, start_response)


app = configure_flask_app(create_flask_app())
app.wsgi_app = InterceptRequestMiddleware(app.wsgi_app)


@app.teardown_appcontext
def close_db(e=None):
    g.pop("db_session", None)


with open(os.path.join(os.path.dirname(__file__), "index.html")) as swagger_ui_file_object:
    swagger_ui_html = swagger_ui_file_object.read()


@app.route("/", methods=["GET", "HEAD"])
def serve_swagger_ui():
    return Response(swagger_ui_html, mimetype="text/html")


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
