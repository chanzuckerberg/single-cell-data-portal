#!/usr/bin/env python

"""
Entry point for starting a local test API server.
"""

import argparse
import logging
import os
import sys

from chalice.cli import CLIFactory

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--host", default="localhost")
parser.add_argument("--port", type=int, default=5000)
parser.add_argument("--no-debug", dest="debug", action="store_false", help="Disable Chalice/Connexion/Flask debug mode")
parser.add_argument("--project-dir", help=argparse.SUPPRESS, default=os.path.join(os.path.dirname(__file__)))
parser.add_argument(
    "--log-level",
    help=str([logging.getLevelName(i) for i in range(0, 60, 10)]),
    choices={logging.getLevelName(i) for i in range(0, 60, 10)},
    default=logging.DEBUG,
)
args = parser.parse_args()
logging.basicConfig(level=args.log_level, stream=sys.stderr)
logging.getLogger("botocore").setLevel(logging.INFO)
factory = CLIFactory(project_dir=args.project_dir, debug=args.debug)

os.environ["DEPLOYMENT_STAGE"] = os.getenv("DEPLOYMENT_STAGE", "test")

# The following corpora snippet is basically stolen from chalice/__init__py:run_local_server
config = factory.create_config_obj(chalice_stage_name=os.environ["DEPLOYMENT_STAGE"])
app_obj = factory.load_chalice_app()
server = factory.create_local_server(app_obj=app_obj, config=config, host=args.host, port=args.port)
server.serve_forever()
