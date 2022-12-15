FROM ubuntu:22.04

ENV APP_NAME=corpora-api
ENV DEPLOYMENT_STAGE=test
ENV EXPORT_ENV_VARS_TO_LAMBDA="APP_NAME DEPLOYMENT_STAGE"
ENV LC_ALL=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive
ENV DD_GEVENT_PATCH_ALL=true

RUN apt-get update && \
    apt-get install -y python3 libhdf5-dev python3-h5py gettext moreutils build-essential libxml2-dev python3-dev python3-pip zlib1g-dev python3-requests python3-aiohttp llvm jq git && \
    rm -rf /var/lib/apt/lists/*

# Don't re-run pip install unless either requirements.txt has changed.
WORKDIR /single-cell-data-portal
ADD requirements.txt requirements-base.txt
ADD backend/api_server/requirements.txt requirements-api.txt


# ddtrace is for Datadog APM metric reporting
RUN python3 -m pip install cmake
RUN python3 -m pip install -r requirements-base.txt -r requirements-api.txt ddtrace
EXPOSE 5000

# Install utilities to /single-cell-data-portal so we can run db migrations.
ADD tests /single-cell-data-portal/tests
ADD scripts /single-cell-data-portal/scripts
ADD backend /single-cell-data-portal/backend

ARG HAPPY_BRANCH="unknown"
ARG HAPPY_COMMIT=""
LABEL branch=${HAPPY_BRANCH}
LABEL commit=${HAPPY_COMMIT}
ENV COMMIT_SHA=${HAPPY_COMMIT}
ENV COMMIT_BRANCH=${HAPPY_BRANCH}

# For Datadog <-> gunicorn integration
# https://docs.datadoghq.com/containers/docker/integrations/?tab=docker#configuration
# https://docs.datadoghq.com/integrations/gunicorn/#metric-collection
LABEL "com.datadoghq.ad.check_names"='["gunicorn"]'
LABEL "com.datadoghq.ad.init_configs"='[{}]'
LABEL "com.datadoghq.ad.instances"='[{ "proc_name": "backend.api_server.app:app" }]'


ADD entrypoint.sh entrypoint.sh
RUN chmod +x entrypoint.sh
CMD ["/single-cell-data-portal/entrypoint.sh"]

# Note: Using just 1 worker for dev/test env. Multiple workers are used in deployment envs, as defined in Terraform code.
# gunicorn --worker-class gevent --workers 1 --bind 0.0.0.0:5000 backend.api_server.app:app --max-requests 10000 --timeout 180 --keep-alive 5 --log-level info
