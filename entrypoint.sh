#!/bin/sh

echo "Starting ddtrace and gunicorn"
export DD_GEVENT_PATCH_ALL=1
ddtrace-run gunicorn --worker-class=gthread --statsd-host=localhost:8125 --workers 1 --bind 0.0.0.0:5000 backend.api_server.app:app --max-requests 10000 --timeout 180 --keep-alive 5 --log-level info --reload