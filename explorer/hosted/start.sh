#!/bin/bash

nginx

gunicorn --worker-class gevent --bind 0.0.0.0:4555 server.ecs.app:application --timeout 60