FROM ubuntu:22.04

ENV APP_NAME=corpora-api
ENV DEPLOYMENT_STAGE=test
ENV LC_ALL=C.UTF-8

RUN apt-get update && \
    apt-get install -y python3 libhdf5-dev python3-h5py gettext moreutils build-essential libxml2-dev python3-dev \
    python3-pip zlib1g-dev python3-requests python3-aiohttp llvm jq && \
    rm -rf /var/lib/apt/lists/*

ADD mock_auth0_pkce.py mock_auth0_pkce.py
ADD requirements.txt requirements.txt
RUN python3 -m pip install -r requirements.txt
EXPOSE 443

ADD container_init.sh container_init.sh

ENTRYPOINT ["./container_init.sh"]
