# Chalice dockerfile
FROM ubuntu:focal-20210119

ENV APP_NAME=corpora-api
ENV DEPLOYMENT_STAGE=dev
ENV EXPORT_ENV_VARS_TO_LAMBDA="APP_NAME DEPLOYMENT_STAGE"
ENV LC_ALL=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y gettext moreutils build-essential libxml2-dev python3-dev python3-pip zlib1g-dev python3-requests python3-aiohttp llvm jq && \
    rm -rf /var/lib/apt/lists/*

# Make python3 the default 'python' executable.
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.8 1

# Don't re-run pip install unless either requirements.txt has changed.
WORKDIR /corpora-data-portal
ADD requirements.txt /corpora-data-portal/requirements.txt
ADD backend/chalice/api_server/requirements.txt /corpora-data-portal/requirements-api.txt
RUN grep -v requirements.txt requirements.txt > reqs.txt \
    && cat requirements-api.txt >> reqs.txt \
    && python3 -m pip install -r reqs.txt

# Install utilities to /corpora-data-portal so we can run db migrations.
ADD tests /corpora-data-portal/tests
ADD scripts /corpora-data-portal/scripts
ADD backend /corpora-data-portal/backend

# Add our api server code as a chalice package in /chalice.
# NOTE: we're relying on .dockerignore to exclude some files
WORKDIR /chalice
ADD backend/chalice/api_server /chalice
RUN mv .chalice/config.json.dev .chalice/config.json

RUN mkdir -p chalicelib config vendor
ADD backend/corpora chalicelib/corpora
ADD backend/config/corpora-api.yml chalicelib/config/corpora-api.yml

CMD ["python3", "/chalice/run_local_server.py"]
