import base64
import hashlib
import logging
import os
import sys
from logging.config import dictConfig
from urllib.parse import urlparse

from flask import json
from flask_cors import CORS
from flask_talisman import Talisman

SERVERDIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(SERVERDIR)


dictConfig(
    {
        "version": 1,
        "formatters": {
            "default": {
                "format": "[%(asctime)s] %(levelname)s in %(module)s: %(message)s",
            }
        },
        "handlers": {
            "wsgi": {
                "class": "logging.StreamHandler",
                "stream": "ext://flask.logging.wsgi_errors_stream",
                "formatter": "default",
            }
        },
        "root": {"level": "INFO", "handlers": ["wsgi"]},
    }
)

try:
    from server.app.app import Server
    from server.common.config.app_config import AppConfig
    from server.common.utils.data_locator import DataLocator, discover_s3_region_name
except Exception:
    logging.critical("Exception importing server modules", exc_info=True)
    sys.exit(1)


class WSGIServer(Server):
    def __init__(self, app_config):
        super().__init__(app_config)

    @staticmethod
    def _before_adding_routes(app, app_config):
        script_hashes = WSGIServer.get_csp_hashes(app, app_config)

        # add the api_base_url to the connect_src csp header.
        extra_connect_src = []
        api_base_url = app_config.server__app__api_base_url
        if api_base_url:
            parse_api_base_url = urlparse(api_base_url)
            extra_connect_src = [f"{parse_api_base_url.scheme}://{parse_api_base_url.netloc}"]

        PLAUSIBLE_URL = "https://plausible.io"

        HUBSPOT_JS_URL = "https://js.hsforms.net"

        HUBSPOT_FORMS_URL = "https://forms.hsforms.com"

        csp = {
            "default-src": ["'self'", HUBSPOT_FORMS_URL, HUBSPOT_JS_URL],
            "form-action": ["'self'", HUBSPOT_FORMS_URL],
            "connect-src": ["'self'", PLAUSIBLE_URL, HUBSPOT_FORMS_URL] + extra_connect_src,
            "script-src": ["'self'", "'unsafe-eval'", PLAUSIBLE_URL, HUBSPOT_FORMS_URL, HUBSPOT_JS_URL] + script_hashes,
            "style-src": ["'self'", "'unsafe-inline'"],
            "img-src": ["'self'", "https://cellxgene.cziscience.com", "data:", HUBSPOT_FORMS_URL],
            "object-src": ["'none'"],
            "base-uri": ["'none'"],
            "frame-ancestors": ["'none'"],
        }

        if not app.debug:
            csp["upgrade-insecure-requests"] = ""

        if app_config.server__app__csp_directives:
            for k, v in app_config.server__app__csp_directives.items():
                csp[k] = csp.get(k, []) + v

        # Add the web_base_url to the CORS header
        web_base_url = app_config.server__app__web_base_url
        if web_base_url:
            web_base_url_parse = urlparse(web_base_url)
            allowed_origins = [f"{web_base_url_parse.scheme}://{web_base_url_parse.netloc}"]
            if os.getenv("DEPLOYMENT_STAGE") in ["Staging", "staging"]:
                allowed_origins.extend(
                    [
                        "https://canary-cellxgene.dev.single-cell.czi.technology/",
                        r"^http://localhost:\d+",
                    ]
                )
            CORS(app, supports_credentials=True, origins=allowed_origins)

        Talisman(
            app,
            force_https=app_config.server__app__force_https,
            frame_options="DENY",
            content_security_policy=csp,
        )

    @staticmethod
    def load_static_csp_hashes(app):
        csp_hashes = None
        try:
            with app.open_resource("../common/web/csp-hashes.json") as f:
                csp_hashes = json.load(f)
        except FileNotFoundError:
            pass
        if not isinstance(csp_hashes, dict):
            csp_hashes = {}
        script_hashes = [f"'{hash}'" for hash in csp_hashes.get("script-hashes", [])]
        if len(script_hashes) == 0:
            logging.error("Content security policy hashes are missing, falling back to unsafe-inline policy")

        return script_hashes

    @staticmethod
    def compute_inline_csp_hashes(app, app_config):
        dataset_configs = [app_config.default_dataset] + list(app_config.dataroot_config.values())
        hashes = []
        for dataset_config in dataset_configs:
            inline_scripts = dataset_config["app"]["inline_scripts"]
            for script in inline_scripts:
                with app.open_resource(f"../common/web/templates/{script}") as f:
                    content = f.read()
                    # we use jinja2 template include, which trims final newline if present.
                    if content[-1] == 0x0A:
                        content = content[0:-1]
                    hash = base64.b64encode(hashlib.sha256(content).digest())
                    hashes.append(f"'sha256-{hash.decode('utf-8')}'")
        return hashes

    @staticmethod
    def get_csp_hashes(app, app_config):
        script_hashes = WSGIServer.load_static_csp_hashes(app)
        script_hashes += WSGIServer.compute_inline_csp_hashes(app, app_config)
        return script_hashes


try:
    app_config = False
    # config file: look first for "config.yaml" in the current working directory
    config_file = "config.yaml"
    config_location = DataLocator(config_file)
    if config_location.exists():
        with config_location.local_handle() as lh:
            logging.info(f"Configuration from {config_file}")
            app_config = AppConfig(lh)
    else:
        # config file: second, use the CXG_CONFIG_FILE
        config_file = os.getenv("CXG_CONFIG_FILE")
        if config_file:
            region_name = discover_s3_region_name(config_file)
            config_location = DataLocator(config_file, region_name)
            if config_location.exists():
                with config_location.local_handle() as lh:
                    logging.info(f"Configuration from {config_file}")
                    app_config = AppConfig(lh)
            else:
                logging.critical(f"Configuration file not found {config_file}")
                sys.exit(1)

    if not app_config:
        logging.critical("No config file found")
        sys.exit(1)

    dataroot = os.getenv("CXG_DATAROOT")
    if dataroot:
        logging.info("Configuration from CXG_DATAROOT")
        app_config.update_server_config(multi_dataset__dataroot=dataroot)

    # complete config
    app_config.complete_config()

    server = WSGIServer(app_config)
    debug = False
    application = server.app

except Exception:
    logging.critical("Caught exception during initialization", exc_info=True)
    sys.exit(1)
    logging.info(f"starting server with multi_dataset__dataroot={app_config.server__multi_dataset__dataroot}")

if __name__ == "__main__":
    try:
        application.run(host=app_config.server__app__host, debug=debug, threaded=not debug, use_debugger=False)
    except Exception:
        logging.critical("Caught exception during initialization", exc_info=True)
        sys.exit(1)
