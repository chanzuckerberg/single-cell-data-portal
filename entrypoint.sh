#!/bin/sh

export DD_AGENT_HOST=$(curl http://169.254.169.254/latest/meta-data/local-ipv4); ddtrace-run python single-cell-data-portal

gunicorn --worker-class gevent --workers 1 --bind 0.0.0.0:5000 backend.api_server.app:app --max-requests 10000 --timeout 180 --keep-alive 5 --log-level info
