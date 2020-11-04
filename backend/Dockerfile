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

# Don't re-run pip install unless requirements.txt has changed.
WORKDIR /chalice
ADD chalice/api_server/requirements.txt /chalice/requirements.txt
RUN python3 -m pip install -r requirements.txt

# Add our api server code.
# NOTE: we're relying on .dockerignore to exclude some files
ADD chalice/api_server /chalice
RUN mv .chalice/config.json.dev .chalice/config.json

RUN mkdir -p chalicelib config vendor
ADD corpora chalicelib/corpora
ADD config/corpora-api.yml chalicelib/config/corpora-api.yml
# Make python3 the default 'python' executable.
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.8 1

CMD ["python3", "run_local_server.py"]
