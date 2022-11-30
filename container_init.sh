#!/bin/bash
if [ "${DEPLOYMENT_STAGE}" == "test" ]; then
  # Use locally-generated cert for HTTPS in containerized local environment; deployed envs use ELB
  echo "--certfile /tmp/pkcs12/server.crt --keyfile /tmp/pkcs12/server.key"
fi
