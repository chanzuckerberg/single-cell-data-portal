#!/usr/bin/env bash
set -e

client_port=${CXG_CLIENT_PORT:-`jq -r '.CXG_CLIENT_PORT' environment.default.json`}
echo -n "localhost:${client_port}/" > .test_base_url.txt

PROJECT_ROOT=$(git rev-parse --show-toplevel)

export PYTHONPATH=${PROJECT_ROOT}  # permits module discovery when run from somewhere other than top level dir

PROJECT_ROOT=${PROJECT_ROOT} python -m server.cli.launch -d $@

