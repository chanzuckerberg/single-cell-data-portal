from collections import defaultdict

import connexion
import json
import logging
import os
import re
import sys
from connexion import FlaskApi, ProblemException, problem
from flask import g, Response
from flask_cors import CORS
from urllib.parse import urlparse

import backend
from backend.corpora.common.utils.json import CustomJSONEncoder
from backend.corpora.common.utils.aws import AwsSecret
from backend.corpora.common.utils.db_session import db_session_manager



def create_flask_app():
    app = connexion.FlaskApp(f"{os.environ['APP_NAME']}-{os.environ['DEPLOYMENT_STAGE']}")
    swagger_spec_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(backend.__file__))), "config", f"{os.environ['APP_NAME']}.yml")
    app.add_api(swagger_spec_path, validate_responses=True)
    return app.app


def configure_flask_app(flask_app):
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
    flask_app.debug = True
    flask_app.logger.setLevel(logging.DEBUG)
    return flask_app


class DatabaseMiddleware:
    def __init__(self, app):
        self.app = app

    def __call__(self, *args, **kwargs):
        g.db_session = db_session_manager()
        g.db_session.__enter__()
        return self.app(*args, **kwargs)


app = configure_flask_app(create_flask_app())
app.wsgi_app = DatabaseMiddleware(app.wsgi_app)


@app.teardown_appcontext
def close_db(e=None):
    g.pop("db_session", None)


@app.teardown_request
def close_transaction(e=None):
    g.db_session.__exit__(*sys.exc_info())


with open(os.path.join(os.path.dirname(__file__), "index.html")) as swagger_ui_file_object:
    swagger_ui_html = swagger_ui_file_object.read()


@app.route("/", methods=["GET", "HEAD"])
def serve_swagger_ui():
    return flask.Response(swagger_ui_html, mimetype="text/html")


@app.errorhandler(ProblemException)
def handle_corpora_error(exception):
    return FlaskApi.get_response(problem(
        exception.status,
        exception.title,
        exception.detail,
        exception.type,
        exception.instance,
        exception.headers,
        exception.ext,
    ))


if __name__ == "__main__":
   app.run(host='0.0.0.0', debug=True)
