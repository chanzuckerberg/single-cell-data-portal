#!/usr/bin/env python3

"""
Entry point for starting a local test API server.
"""

import argparse
import functools
import logging
import os
import sys
import time

from chalice.deploy.validate import validate_routes
from chalice.cli import CLIFactory, reloader


def get_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="localhost")
    parser.add_argument("--port", type=int, default=5000)
    parser.add_argument(
        "--no-debug", dest="debug", action="store_false", help="Disable Chalice/Connexion/Flask debug mode"
    )
    parser.add_argument("--project-dir", help=argparse.SUPPRESS, default=os.path.join(os.path.dirname(__file__)))
    parser.add_argument(
        "--log-level",
        help=str([logging.getLevelName(i) for i in range(0, 60, 10)]),
        choices={logging.getLevelName(i) for i in range(0, 60, 10)},
        default=logging.INFO,
    )
    os.environ["DEPLOYMENT_STAGE"] = os.getenv("DEPLOYMENT_STAGE", "test")
    stage = os.environ.get("DEPLOYMENT_STAGE")

    args = parser.parse_args()
    return args, stage


def run_server(args, stage):
    logging.basicConfig(level=args.log_level, stream=sys.stderr)
    logging.getLogger("botocore").setLevel(logging.INFO)
    factory = CLIFactory(project_dir=args.project_dir, debug=args.debug)

    # The following code snippet is basically stolen from chalice/__init__py:local
    config = factory.create_config_obj(chalice_stage_name=stage)
    app_obj = factory.load_chalice_app()

    # We don't create the server here because that will bind the
    # socket and we only want to do this in the worker process.
    server_factory = functools.partial(create_local_server, factory, config, app_obj, args.host, args.port, stage)

    # support autoreload
    project_dir = config.project_dir
    print(sys.argv)
    print(os.environ)
    rc = reloader.run_with_reloader(server_factory, os.environ, project_dir)
    # Click doesn't sys.exit() with the RC this function.  The
    # recommended way to do this is to use sys.exit() directly,
    # see: https://github.com/pallets/click/issues/747
    sys.exit(rc)


def create_local_server(factory, config, app, host, port, stage):
    # Check that `chalice deploy` would let us deploy these routes, otherwise
    # there is no point in testing locally.
    routes = config.chalice_app.routes
    validate_routes(routes)
    server = factory.create_local_server(app, config, host, port)
    return server


if __name__ == "__main__":
    args, stage = get_args()
    if os.getenv("RESTART_ON_FAILURE"):
        while True:
            try:
                run_server(args, stage)
            except Exception as err:
                logging.exception(err)
                time.sleep(1)
    else:
        run_server(args, stage)
