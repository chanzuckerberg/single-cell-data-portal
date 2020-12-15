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
import signal

from chalice.deploy.validate import validate_routes
from chalice.cli import CLIFactory, reloader
from chalice.local import LocalDevServer

from six.moves.socketserver import ForkingMixIn
from six.moves.BaseHTTPServer import HTTPServer


class ForkedHTTPServer(ForkingMixIn, HTTPServer):
    """Forking mixin to better support browsers.

    When a browser sends a GET request to Chalice it keeps the connection open
    for reuse. In the single threaded model this causes Chalice local to become
    unresponsive to all clients other than that browser socket. Even sending a
    header requesting that the client close the connection is not good enough,
    the browswer will simply open another one and sit on it.
    """
    allow_reuse_address = True
    timeout = 2

    # Make sure to reap child processes.
    def shutdown(self):
        if self.active_children is None:
            return
        for child_pid in self.active_children:
            os.kill(child_pid, signal.SIGTERM)
        self.server_close()


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

    # handle sigterm
    signal.signal(signal.SIGTERM, lambda *args: sys.exit(0))

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
    server = LocalDevServer(app, config, host, port, server_cls=ForkedHTTPServer)
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
