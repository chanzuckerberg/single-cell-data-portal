import json
import os
import threading
import time
from urllib.parse import urlparse

from connexion import FlaskApi, FlaskApp, ProblemException, problem
from flask import Response, g, request
from flask_cors import CORS
from memory_profiler import memory_usage
from server_timing import Timing as ServerTiming
from swagger_ui_bundle import swagger_ui_path

from backend.common.server.datadog import initialize_datadog_tracing
from backend.common.server.logger import configure_logging
from backend.common.server.request_id import generate_request_id, get_request_id
from backend.common.utils.aws import AwsSecret
from backend.common.utils.json import CustomJSONEncoder

DEPLOYMENT_STAGE = os.environ["DEPLOYMENT_STAGE"]
APP_NAME = "{}-{}".format(os.environ.get("APP_NAME", "api"), DEPLOYMENT_STAGE)


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


def register_routes(app, api_base_paths):
    def apis_landing_page() -> str:
        """
        Render a page that displays links to all APIs
        """
        links = [f'<a href="{base_path}/ui/">{base_path}</a></br>' for base_path in api_base_paths]
        return f"""
        <html>
        <head><title>cellxgene Platform APIs</title></head>
        <body><h1>cellxgene Platform APIs</h1>{"".join(links)}</body>
        </html>
        """

    def before_request():
        g.start_time = time.time()
        g.request_id = generate_request_id()
        g.track_memory = request.headers.get("X-Enable-Memory-Tracking", "false").lower() == "true"
        if g.track_memory:
            g.memory_usage = []
            g.memory_thread = threading.Thread(target=_track_memory_usage, args=(g,))
            g.memory_thread.start()
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

    def after_request(response: Response):
        if g.track_memory:
            g.track_memory = False
            g.memory_thread.join()
            peak_memory = max(g.memory_usage) if g.memory_usage else 0
            response.headers.add("X-Peak-Memory-Usage", f"{peak_memory:.3f} MiB")

        total_time = time.time() - g.start_time
        response.headers.add("X-Request-Id", get_request_id())
        response.headers.add("X-Total-Request-Time", f"{total_time:.3f}s")

        app.logger.info(
            dict(
                type="RESPONSE",
                details=dict(
                    status_code=response.status_code,
                    content_length=response.content_length,
                    response_time=total_time,
                ),
            )
        )
        return response

    def _track_memory_usage(g_context):
        while g_context.track_memory:
            usage = memory_usage(-1, interval=0.25, timeout=0.25)
            if usage:
                g_context.memory_usage.append(usage[0])

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

    def handle_internal_server_error(exception):
        app.logger.exception("InternalServerError", exc_info=exception)
        return FlaskApi.get_response(problem(500, "Internal Server Error", "Internal Server Error"))

    app.add_url_rule("/", "apis_landing_page", apis_landing_page)
    app.before_request(before_request)
    app.after_request(after_request)
    app.register_error_handler(ProblemException, handle_corpora_error)
    app.register_error_handler(Exception, handle_internal_server_error)

    if DEPLOYMENT_STAGE == "test":

        def raise_exception():
            raise Exception("testing")

        app.add_url_rule("/exception", "raise_exception", raise_exception)


def create_api_app(api_paths_and_spec_files, **server_args):
    """
    Creates and configures a Connexion API application.

    This function initializes a Connexion application with the given API paths and specification files.
    It sets up the application with necessary configurations, including URL rules, request handlers,
    and error handlers. The function also configures the Flask application with appropriate settings
    based on the deployment stage.

    Args:
        api_paths_and_spec_files (list of tuples): A list where each tuple contains a base path and a specification file.
        **server_args: Additional arguments to be passed to the FlaskApp.

    Returns:
        flask.Flask: The configured Flask application wrapped by Connexion.
    """

    connexion_app = FlaskApp(APP_NAME, specification_dir="backend", server_args=server_args)

    # From https://github.com/zalando/connexion/issues/346
    connexion_app.app.url_map.strict_slashes = False

    def add_api(base_path, spec_file):
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

    apis = {}
    for base_path, spec_file in api_paths_and_spec_files:
        apis[base_path] = add_api(base_path=base_path, spec_file=spec_file)

    app = connexion_app.app

    # store the apis in the connexion app
    app.apis = apis

    app = configure_flask_app(app)
    register_routes(app, [base_path for base_path, _ in api_paths_and_spec_files])
    return app


configure_logging(APP_NAME)
initialize_datadog_tracing()
