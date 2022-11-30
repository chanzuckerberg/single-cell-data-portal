#!/bin/bash
echo
echo " ====="
echo "| starting backend container"
echo " ====="
echo

if [ "${DEPLOYMENT_STAGE}" == "test" ]; then
  # Use locally-generated cert for HTTPS in containerized local environment; deployed envs use ELB
  echo
  echo " ====="
  echo "| Using test cert"
  echo " ====="
  echo
  export HTTPS_CERT_AND_KEY="--certfile /tmp/pkcs12/server.crt --keyfile /tmp/pkcs12/server.key"
fi
