#!/bin/sh

echo "Starting ddtrace and gunicorn"
export DD_GEVENT_PATCH_ALL=1
# --statsd-host and --name are for https://docs.datadoghq.com/integrations/gunicorn/#metric-collection
# --name must agree with Dockerfile's `LABEL "com.datadoghq.ad.instances"` line
ddtrace-run gunicorn --workers 1 --bind 0.0.0.0:5000 backend.api_server.app:app --max-requests 10000 --timeout 180 --keep-alive 5 --log-level info --reload --worker-class gevent --preload --statsd-host localhost:8125 --name backend
