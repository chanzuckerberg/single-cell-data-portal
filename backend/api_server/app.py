import json
import os
import time
from urllib.parse import urlparse

import connexion
from connexion import FlaskApi, ProblemException, problem
from ddtrace import patch_all, tracer
from ddtrace.filters import FilterRequestsOnUrl
from flask import Response, g, request
from flask_cors import CORS
from server_timing import Timing as ServerTiming
from swagger_ui_bundle import swagger_ui_path

from backend.api_server.logger import configure_logging
from backend.api_server.request_id import generate_request_id, get_request_id
from backend.common.utils.aws import AwsSecret
from backend.common.utils.json import CurationJSONEncoder, CustomJSONEncoder
from backend.gene_info.api.ensembl_ids import GeneChecker

DEPLOYMENT_STAGE = os.environ["DEPLOYMENT_STAGE"]
APP_NAME = "{}-{}".format(os.environ.get("APP_NAME", "api"), DEPLOYMENT_STAGE)


configure_logging(APP_NAME)


def should_configure_datadog_tracing():
    return (
        DEPLOYMENT_STAGE in ["dev", "staging", "prod"]
        and os.environ.get("DD_AGENT_HOST", None)
        and os.environ.get("DD_TRACE_AGENT_PORT", None)
    )


if should_configure_datadog_tracing():
    # Datadog APM tracing
    # See https://ddtrace.readthedocs.io/en/stable/basic_usage.html#patch-all

    tracer.configure(
        hostname=os.environ["DD_AGENT_HOST"],
        port=os.environ["DD_TRACE_AGENT_PORT"],
        # Filter out health check endpoint (index page: '/')
        settings={
            "FILTERS": [
                FilterRequestsOnUrl([r"http://.*/$"]),
            ],
        },
    )
    patch_all()


def create_flask_app():
    connexion_app = connexion.FlaskApp(
        APP_NAME, specification_dir="backend", server_args=dict(static_folder="backend/api_server/static")
    )

    # From https://github.com/zalando/connexion/issues/346
    connexion_app.app.url_map.strict_slashes = False

    def add_api(base_path, spec_file):
        api_base_paths.append(base_path)
        api = connexion_app.add_api(
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
        return api

    add_api(base_path="/dp", spec_file="portal/api/portal-api.yml")
    curation_api = add_api(base_path="/curation", spec_file="curation/api/curation-api.yml")
    curation_api.blueprint.json_encoder = CurationJSONEncoder
    add_api(base_path="/wmg", spec_file="wmg/api/wmg-api.yml")
    add_api(base_path="/wmg/v2", spec_file="wmg/api/wmg-api-v2.yml")
    add_api(base_path="/gene_info", spec_file="gene_info/api/gene-info-api.yml")

    # Initialize gene checker to go ahead and create a dictionary of all
    # gene ensembl ID to gene name mappings for the gene_info API endpoint.
    # This is done here because the initialization of gene checker takes some
    # time and, if done on API request, would significantly slow performance.
    GeneChecker()

    return connexion_app.app


def configure_flask_app(flask_app):
    flask_app.debug = DEPLOYMENT_STAGE != "prod"

    # set the flask secret key, needed for session cookies
    flask_secret_key = "OpenSesame"
    allowed_origins = []
    deployment_stage = os.environ["DEPLOYMENT_STAGE"]
    if deployment_stage not in ["prod"]:
        allowed_origins.extend([r"https?://.*\.corporanet\.local:\d+", r"^https?://localhost:\d+"])
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
    ServerTiming(flask_app, force_debug=True)
    # FIXME, enforce that the flask_secret_key is found once all secrets are setup for all environments
    require_secure_cookies = not bool(os.getenv("DEV_MODE_COOKIES"))
    flask_app.config.update(
        SECRET_KEY=flask_secret_key,
        SESSION_COOKIE_SECURE=require_secure_cookies,
        SESSION_COOKIE_HTTPONLY=True,
        SESSION_COOKIE_SAMESITE="Lax",
        JSON_SORT_KEYS=True,
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


if DEPLOYMENT_STAGE == "test":

    @app.route("/exception")
    def raise_exception():
        raise Exception("testing")


@app.before_request
def before_request():
    g.start = time.time()
    g.request_id = generate_request_id()
    app.logger.info(
        dict(
            type="REQUEST",
            details=dict(
                url=request.path,
                method=request.method,
                content_length=request.content_length,
            ),
        )
    )


@app.after_request
def after_request(response: Response):
    app.logger.info(
        dict(
            type="RESPONSE",
            details=dict(
                status_code=response.status_code,
                content_length=response.content_length,
                response_time=time.time() - g.start,
            ),
        )
    )
    response.headers["X-Request-Id"] = get_request_id()
    return response


@app.errorhandler(ProblemException)
def handle_corpora_error(exception):
    if exception.status >= 500:
        app.logger.error("InternalServerError", exc_info=exception)
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


@app.errorhandler(Exception)
def handle_internal_server_error(exception):
    app.logger.exception("InternalServerError", exc_info=exception)
    return FlaskApi.get_response(problem(500, "Internal Server Error", "Internal Server Error"))


if __name__ == "__main__":
    app.run(host="0.0.0.0", debug=True)
