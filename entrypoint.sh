#!/bin/sh

echo "Starting ddtrace"
ddtrace-run python single-cell-data-portal

echo "Starting gunicorn"
gunicorn --worker-class gevent --workers 1 --bind 0.0.0.0:5000 backend.api_server.app:app --max-requests 10000 --timeout 180 --keep-alive 5 --log-level info
