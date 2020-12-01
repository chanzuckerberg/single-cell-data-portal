# Chalice dockerfile
FROM ubuntu:focal

ENV APP_NAME=corpora-api
ENV DEPLOYMENT_STAGE=dev
ENV EXPORT_ENV_VARS_TO_LAMBDA="APP_NAME DEPLOYMENT_STAGE"
ENV LC_ALL=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y gettext moreutils build-essential libxml2-dev python3-dev python3-pip zlib1g-dev python3-requests python3-aiohttp llvm

# RUN envsubst < config/iam-policy-templates/corpora-api-lambda.json > .chalice/policy-$(DEPLOYMENT_STAGE).json

ADD tests /tests
ADD scripts /scripts
WORKDIR /backend
ADD backend /backend
ADD requirements.txt /requirements.txt

# Don't re-run pip install unless requirements.txt has changed.
RUN python3 -m pip install -r /requirements.txt

RUN mv chalice/api_server/.chalice/config.json.dev chalice/api_server/.chalice/config.json

RUN mkdir -p chalicelib config vendor
ADD backend/corpora chalicelib/corpora
ADD backend/config/corpora-api.yml chalicelib/config/corpora-api.yml
# Make python3 the default 'python' executable.
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.8 1

CMD ["python3", "/backend/chalice/api_server/run_local_server.py"]
