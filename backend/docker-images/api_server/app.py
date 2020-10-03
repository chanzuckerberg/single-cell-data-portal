import json
import os
import re
from collections import defaultdict

import connexion
from corpora.common.utils.aws_secret import AwsSecret
from corpora.common.utils.json import CustomJSONEncoder
from flask import make_response
from flask_cors import CORS


def remove_body_from_api_request_or_response_for_logging(logging_message):
    return logging_message.to_dict().pop("body", None)


class CorporaServer():
    def __init__(self):
        self.deployment_stage = os.environ["DEPLOYMENT_STAGE"]

        flask_api = connexion.FlaskApp(__name__, specification_dir='config/')
        flask_api.add_api(f"{os.environ['APP_NAME']}.yml", validate_responses=True)
        self.corpora_app = flask_api.app

        # Allow CORS for all domains on all routes
        CORS(self.corpora_app)

        self.configure_server_app()

    def configure_server_app(self):
        self.corpora_app.json_encoder = CustomJSONEncoder

        self.setup_database_connection_for_app()
        self.setup_authentication_for_app()
        self.setup_api_resources()

    def setup_application_logging(self):
        if self.deployment_stage != "prod":
            self.corpora_app.debug = True

    def setup_database_connection_for_app(self):
        database_secret_name = f"corpora/backend/{self.deployment_stage}/database"
        database_location = json.loads(AwsSecret(database_secret_name).value)
        self.corpora_app.config['SQLALCHEMY_DATABASE_URI'] = database_location

    def setup_authentication_for_app(self):
        auth0_secret_name = f"corpora/backend/{self.deployment_stage}/auth0-secret"
        auth0_secret = json.loads(AwsSecret(auth0_secret_name).value)
        if auth0_secret:
            flask_secret_key = auth0_secret.get("flask_secret_key", "OpenSesame")  # Arbitrary secret if none is set.
        else:
            flask_secret_key = "OpenSesame"
        self.corpora_app.config['SECRET_KEY'] = flask_secret_key

    def setup_api_resources(self):
        def dispatch(*args, **kwargs):
            self.corpora_app.log.info(
                f"Request: {remove_body_from_api_request_or_response_for_logging(self.corpora_app.current_request)}")

            uri_params = self.corpora_app.current_request.uri_params or {}
            resource_path = self.corpora_app.current_request.context["resourcePath"].format(**uri_params)
            req_body = self.corpora_app.current_request.raw_body if self.corpora_app.current_request._body is not \
                                                                    None else None

            with self.corpora_app.test_request_context(
                    path=resource_path,
                    base_url="https://{}".format(self.corpora_app.current_request.headers["host"]),
                    query_string=self.corpora_app.current_request.query_params,
                    method=self.corpora_app.current_request.method,
                    headers=list(self.corpora_app.current_request.headers.items()),
                    data=req_body,
                    environ_base=self.corpora_app.current_request.stage_vars,
            ):
                flask_response = self.corpora_app.full_dispatch_request()

            self.corpora_app.current_request.query_params.log.info(
                f"Response: {remove_body_from_api_request_or_response_for_logging(flask_response)}")

            return flask_response

        routes = defaultdict(list)
        for rule in self.corpora_app.url_map.iter_rules():
            routes[re.sub(r"<(.+?)(:.+?)?>", r"{\1}", rule.rule).rstrip("/")] += rule.methods
        for route, methods in routes.items():
            self.corpora_app.route(route, methods=list(set(methods) - {"OPTIONS"}))(dispatch)

        # Add top level route to swagger page
        self.corpora_app.add_url_rule("/", "swagger", self.serve_swagger_ui)

    def serve_swagger_ui(self):
        with open("index.html") as swagger_ui_file_object:
            return make_response(swagger_ui_file_object.read())


if __name__ == "__main__":
    corpora_server = CorporaServer()
    corpora_server.corpora_app.run(debug=True, host='0.0.0.0')
